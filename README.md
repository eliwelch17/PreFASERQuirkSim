# PreFASERQuirkSim
Simulation of Quirks from ATLAS IP toward FASER adapted from MATHEMATICA simulation written by Junle Pei and Jinmian Li used in https://arxiv.org/pdf/2404.13814 

In an attempt to make this as simple as possible for those who may use this (including theorists without extensive practice with cmake), I have tried to contain the program to mostly a single script with a single compile line. If more is needed, a proper cmake build may be used in the future.



## Requirements: 
    g++ v7 or later (check: g++ --version), c++17

## Compile: 
    g++ -std=c++17 -o quirk_run quirk_run.cxx src/KDTree3D.cxx

## Usage: 
    
    ./quirk_run [-b <back_value m>] [-l <lambda_value eV>] [-s <seed>] [-n <# quirks>] [-d <stepsize divider>] [-t (trajectory output flag)] <input file>


    -b how far back to simulate the quirks in m (defult is 480m)

    -l energy scale lambda in eV

    -s random seed

    -n number of quirks to simualte (default is all)

    -d: the time step size is propotional to lambda^2/d (default is 10000)

    -t: option to record the position of the quirks every .01 ns simluated time, given for quirk1 (default is false)


## Validation: 
    The code validates well against the original Mathematica code. The remaining issue is a small numerical instability resulting from loss of precision in Mathematica which compounds over many steps, which means it becomes noticeable for large lambda. Even for large lambda however, the arrival time and maximum transverse position and momentum remain unchanged. 

    Comparrisons of the quirks trajectories for various lambdas and 200 GeV quirks are located in validation_plots folder. Note that the sampling frequency of the trajectory is less than the oscillation frequency.
    
    NB: It is rather difficult to get Mathematica and C++ to extract random numbers with the same seed so for the validation runs, the Gaussian de/dx functions simply returned a constant to be able to compare the two. 

    Regardless, the differences accrued by numerical instability are a negligible error compared to those introduced by the random walk effect from de/dx.

## Simulation time improvement:
    the c++ code seems to be about 15-80x faster than the Mathematica code

## Magnetic field

    The original mathematica code only included the D1 and D2 LHC dipole magnets and provided a nearly uniform field. This c++ simulation now includes realistic field maps of each magnet: D1,D2, MB, MQAB2 and MQX13. The quadripole magnets are to be added soon. A KDTree is used to quickly obtain the field value at a given point.


    Currently the B-field resolution has 1/5 x 1/5  (transverse plane) of that in the provided field maps.

   






