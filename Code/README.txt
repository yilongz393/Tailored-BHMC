# Global Energy Minimum Optimization of Interfacial Assembly of Binary Nanoparticles (NPs)

* * *
The attached code implements a basin-hopping Monte Carlo (BHMC) algorithm specifically tailored for determining the ground-state configuration of NPs trapped at a fluid-fluid interface. The algorithm begins by randomly generating an initial seed configuration of NPs within a cuboid box in which the interfacial plane is fixed at z=0. One of the six MC moves (based on a sequence) is then implemented to produce a new configuration of NPs, which is further optimized to it nearest local energy minimum by using local minimization approaches. In particular, the Polak-Ribiere variant of the conjugate gradient algorithm is used for local minimization, where an inexact line search is performed based on the Wolfe conditions. The Metropolis criterion is then used to accept or reject this new NP configuration, which is then set as the seed structure for the next iteration if accepted. 
* * * 

## Introduction
* * *
1. The algorithm is implemented in C++, which consists of two components: "main.cpp" and "BHMC.h". The header file "BHMC.h" is a class for the tailored basin-hopping Monte Carlo algorithm, where the system parameters are first defined, the functions for generating initial configuration of NPs and for calculating the inter-particle and interfacial forces are then generated, and the functions associated with conjugate gradient minimization and MC moves are also included. The "main.cpp" performed the operations of this BHMC algorithm. 
2. The code takes an input file "Input", which includes 8 system parameters in total, one per row, as specified below
* * *
first row: Number of type I NPs
second row: Number of type II NPs
third row: Radius of NP I
fourth row: Radius of NP II
fifth row: Omega --- the ratio of the differences in the surface energy of the two types of NPs with the two layers
sixth row: Eff --- factor accounting for the effective radius of the surface-functionalized NPs, which was set to 1.0 throughout the study
seventh row: Chi --- ratio of the difference in the surface energy of NP I with the two fluids to the surface tension between the two fluids
eighth row: Epsilon --- the relative strength of inter-particle interactions
* * *

## Code Execution 
* * *
1. Download "main.cpp", "BHMC.h", "Makefile", and "Input" and navigate to the folder including these files
2. In the terminal, type
*
$ make
*
to generate a compiled executable code named "BHMC". It should take ~10 seconds to compile the code. To run the code, type in the terminal
*
$ ./BHMC Input
*
The code performs roughly 25 MC steps per second for a 12-particle system. Typical run times till convergence can vary from 1 minute to 100 hours, depending on system size and complexity.
* * * 

## Output Files
* * * 
The code generates 4 output files:
1. GM_ENERGIES.DAT --- Current best energy minima (till that MC step), where the last one indicates the best obtained energy minimum optimized by the algorithm
2. GM_MOVIE.xyz --- NP configurations in xyz format corresponding to the energy minima in the file "GM_ENERGIES.DAT", where the last frame indicates the best NP structure optimized by the algorithm
3. LM_ENERGIES.DAT --- All accepted local energy minima
4. LM_MOVIE.xyz --- NP configurations in xyz format corresponding to the accepted local energy minima in the file "LM_ENERGIES.DAT"
* * *
In addition to these 4 files, the code also prints on the screen information about initializing NP configuration, conjugate gradient operation, and Monte Carlo moves
* * *

## Additional Information
* * *
1. The default box dimension is "15*15*10 times (Radius of NP I)". Periodic boundary conditions were applied to x and y dimensions, and fixed boundaries were used in the z dimension. To change these box dimensions (e.g., for larger system size), go to Line 53 in "BHMC.h" and change the parameter initialization of "HEIGHT" (for dimension z) and "DENS" (for dimension x and y)
2. An example for searching the hexagonally arranged bilayer is included
* * *

## If you use this code, please cite our paper: 

Zhou, Y., Arya, G. Discovery of two-dimensional binary nanoparticle superlattices using global Monte Carlo optimization. Nat Commun 13, 7976 (2022). https://doi.org/10.1038/s41467-022-35690-8
