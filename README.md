# PreFASERQuirkSim

PreFASERQuirkSim is a C++ simulation of fermionic, colorless quirk-pair transport from the ATLAS interaction point (IP) to FASER. It is adapted from the Mathematica simulation by Junle Pei and Jinmian Li used in arXiv:2404.13814](https://arxiv.org/abs/2404.13814). Macroscopic quirk oscillations present a unique challenge for long-lived particle simulations. Accurately modeling their transport and interactions with material requires integration steps much shorter than the oscillation period, making simulation runtime a significant bottleneck. This program is a  dedicated C++ simulation to propagate quirks efficiently from the ATLAS interaction point to the FASER detector while preserving the required accuracy in their dynamics and material interactions, and offers signicant simulation time reduction from the initial mathematica simulation.

## Requirements

- A C++ compiler with C++17 support
- `g++` version 7 or later

Check the installed compiler version with:

```bash
g++ --version
```

## Compilation

From the project root, run:

```bash
g++ -std=c++17 -o quirk_run quirk_run.cxx src/*
```

## Usage

```bash
./quirk_run [options] <input-file>
```

### Options

| Option | Description |
| --- | --- |
| `-f <front>` | Starting position in meters. The default is `19`. Transport from the IP at `z = 0` to `front` is handled analytically when the fast-transport condition is met. |
| `-b <back>` | Final longitudinal limit in metres. The default is `474.4`. The simulation returns the quirks at the final minimum in their oscillation before this limit for use by the Athena Geant4 quirks extension. |
| `-l <lambda>` | Confinement scale, Lambda, in eV. The default is `500`. |
| `-betaCut <beta>` | Minimum pair beta below which the event is stopped. The default is `0.1`. |
| `-s <seed>` | Random-number seed. The default is `0`. |
| `-n <count>` | Number of quirk pairs to simulate. By default, all input events are simulated. |
| `-d <divider>` | Time-step divider. The step size is proportional to `Lambda^2 / divider`. The default is `10000`. |
| `-skip <count>` | Number of input events to skip. The default is `0`. |
| `-r <run-number>` | Run number included in the output filename. The default is `0`. |
| `-t` | Write trajectory output sampled every `0.01 ns`. Disabled by default. |

Example:

```bash
./quirk_run -f 19 -b 474.4 -l 500 -s 0 -n 100 Quirk_masses/quirkE_200GeV_0004.dat
```

For Lambda above approximately `1 keV`, numerical instability can become
non-negligible with the default divider. A divider of `15000` to `20000` is
recommended for those runs.

## Material Interactions

During normal step-by-step transport, the simulation determines which material
each quirk occupies at every integration step and applies the corresponding
ionization energy loss. The modeled material regions are:

- TAS copper from approximately `19.0 m` to `20.8 m`, outside its central 17 mm`-radius aperture
- TAN copper from approximately `140.0 m` to `141.0 m`, including its transverse geometry and two `25 mm`-radius beam holes centred at `y = +/-80 mm`
- Concrete from `380 m` to `390 m`
- Rock from `390 m` to `474.6 m`

The deterministic stopping-power calculation is first used to estimate the
distance travelled during the step. The force is then recalculated with a
Gaussian fluctuation in the material energy loss before the particle state is
advanced. Material is handled independently for each member of the quirk pair.

## Fast High-Lambda Transport

For Lambda values of approximately `3 keV` or greater, the study in
[arXiv:2404.13814](https://arxiv.org/abs/2404.13814) finds that deflection from
material and magnetic fields is negligible relative to the strong infracolor
force. These deflections are rapidly washed out over many oscillations, making per-step transport over the entire IP-to-FASER distance unnecessarily expensive.

To use fast transport over the full distance, set `front` equal to `back`. For
example, with the default back position:

```bash
./quirk_run -f 474.4 -b 474.4 -l 3000 <input-file>
```

The simulation analytically transports the pair from `z = 0` to `front` using
its momentum, mass, and Lambda. Although transverse deflections can be neglected in this regime, accumulated ionization loss cannot. The skipped material loss is therefore calculated with precomputed range tables and applied to the quirk momenta.

The range-table correction includes the same copper, concrete, and rock regions
used during normal transport. It calculates an effective path length using the
oscillation factor and TAS/TAN `Loct()` acceptance, obtains the corresponding
beta reduction, and applies a common scale factor to both quirks' three-momenta.
This preserves their momentum sharing while slowing the pair.

## Lifetime and Radiation Estimates

Each output row includes lifetime survival weights and RMS transverse-deflection widths at the first preshower plane (`480.42 m`). Both are calculated for infracolor-gluon radiation probabilities of `epsilon = 0.07, 0.10, 0.13`, giving a central value of `0.10` per oscillation with a `+/-30%` variation. These hard-coded values can be changed in `quirk_run.cxx`; epsilon affects the reported weights and widths but does not alter the transported trajectory. Applying the lifetime as an event weight, rather than modeling decays at each transport step, is a good approximation because the decay-distance distribution is much narrower than the flight distance to FASER. Most quirk pairs therefore either decay well before reaching FASER (`w ~ 0`) or survive with near certainty (`w ~ 1`), energy loss from IC gluon radiaition  is negligible compared to ionization energy loss in the latter case.

## Magnetic Field

The original Mathematica simulation included only the D1 and D2 LHC dipole magnets with a nearly uniform field. This implementation uses realistic field maps for the following magnets:

- D1 and D2 dipoles
- Inner quadrupoles
- Reverse inner quadrupoles
- Main dipoles
- Main quadrupoles

A three-dimensional KD-tree is used to retrieve the magnetic field efficiently at each position. The field maps have `1 cm x 1 cm` resolution.

The inner-quadrupole fields differ between years only by sign, so they are
expected to have the same effect on quirk acceptance at FASER. The 2024 scaling factor is used because that year has the largest integrated luminosity.

Magnetic-field information is courtesy of Alex Keyken and the
[BDSIM development team](https://www.pp.rhul.ac.uk/bdsim/manual/).

## Performance

This C++ implementation has been seen to be approximately 20 to 100 times faster than the original Mathematica simulation, depending on the simulation parameters, quirk model, and quirk trajectory.
