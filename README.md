# PreFASERQuirkSim
Simulation of Quirks from ATLAS IP toward FASER adapted from MATHEMATICA simulation written by Junle Pei and Jinmian Li used in https://arxiv.org/pdf/2404.13814 

## Requirements: 
    g++ v7 or later (check: g++ --version), c++17

## Compile: 
    g++ -std=c++17 -o quirk_run quirk_run.cxx

## Usage: 
    ./quirk_run -l[lambda in eV] -b [final z position in m] [-s seed] <infile>

## Validation: 
    The code currently seems to validate well for the most part against the original Mathematica code. The remaining issue is a small numerical instability resulting from loss of precision in Mathematica which compounds over many steps, which means it becomes noticeable for large lambda. Even for large lambda however, the arrival time and maximum transverse position and momentum remain unchanged. 

    Comparrisons of the quirks trajectories for various lambdas and 200 GeV quirks are located in validation_plots folder. Note that the sampling frequency of the trajectory is less than the oscillation frequency.
    
    NB: It is rather difficult to get Mathematica and C++ to extract random numbers with the same seed so for the validation runs, the Gaussian de/dx functions simply returned a constant to be able to compare the two. 

    Regardless, the differences accrued by numerical instability are a negligible error compared to those introduced by the random walk effect from de/dx.

## Simulation time improvement:
    the c++ code seems to be about 15-80x faster than the Mathematica code

   






