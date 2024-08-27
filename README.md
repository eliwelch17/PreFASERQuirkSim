# PreFASERQuirkSim
Simulation of Quirks from ATLAS IP toward FASER adapted from MATHEMATICA simulation written by Junle Pei and Jinmian Li used in https://arxiv.org/pdf/2404.13814 

Requirements: g++ v7 or later (check: g++ --version)

Compile: g++ -std=c++17 -o quirk_run quirk_run.cxx

Usage: 
    ./quirk_run -l[lambda in eV] -b [final z position in m] [-s seed] <infile>

Validation: 
    The code currently seems to validate well for the most part against the original Mathematica code. The remaining issue is a small numerical instability (I am currently trying to midigate) which compounds over many steps. 
    
    It is rather difficult to get Mathematica and C++ to extract random numbers with the same seed so currently the Gaussian de/dx functions simply return a constant to compare with the Mathematica code.

Simulation time improvement:
    so far the c++ simulation seems to be about 15-80x faster than the Mathematica code

   






