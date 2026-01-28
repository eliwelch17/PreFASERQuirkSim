# Range tables

This directory contains **precomputed range tables** for the quirk ionization-loss model used in `quirk_run.cxx`.

Each file `range_<mass>.h` (mass in GeV) provides:
- `beta_grid[RANGE_N]`: \(\beta\) values from **0.1 → 0.9999**, log-spaced in \((1-\beta)\)
- `range_Cu_um[RANGE_N]`: range in **µm** through copper to slow from \(\beta\) to **0.1**
- `range_Cc_um[RANGE_N]`: same for concrete
- `range_Rock_um[RANGE_N]`: same for rock

## Conventions

- **β cutoff**: tables are constructed with `beta_min = 0.1`. If \(\beta < 0.1\), the simulation should abort the event.
- **Energy-loss kernel**: `EoxCuAll / EoxCcAll / EoxRockAll` return \(S(\beta)\) in **(1e-16 GeV²)** (the internal convention used in `quirk_run.cxx`).
- **Conversion to energy loss over length**:

\[
dE[\mathrm{GeV}] = S(\beta)\cdot dL_{\mu m}\cdot \left(\frac{10^{-16}}{\mathrm{hbar}c}\right)
\]

with \(\mathrm{hbar}c = 1.9732697\times 10^{-10}\,\mathrm{GeV\cdot \mu m}\), so \(\mathrm{conv}=10^{-16}/(\mathrm{hbar}c)=5.0677307\times 10^{-7}\).

## Regenerating

Run:

```bash
python3 PreFASERQuirkSim/tools/generate_range_tables.py
```

Optional arguments:
- `--n 300` (grid size)
- `--mass-min 50 --mass-max 700 --mass-step 25`
- `--beta-min 0.1 --beta-max 0.9999`
- `--outdir PreFASERQuirkSim/rangeTables`


