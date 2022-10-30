# Global Energy Minimum Optimization of Interfacial Assembly of Binary Nanoparticles (NPs)
* * *
The project implemented a basin-hopping Monte Carlo algorithm tailored for interface-trapped structures to rapidly determine the ground-state configuration of NPs. In the algorithm, the initial seed configuration of NPs was randomly generated in a cuboid box in which the interfacial plane was fixed at z=0. A sequence of MC moves was then implemented to produce new seed configurations, which were further optimized by using local minimization approaches. In particular, the Polak-Ribiere variant of the conjugate gradient algorithm was used for local minimization and the inexact line search was performed based on the Wolfe conditions. The Metropolis criterion was then used to update the NP configuration, which was then set as the seed structure for the next iteration. 
* * * 
## Introduction
* * *
1. The algorithm was implemented in C++, which consists of two components: "main.cpp" and "BHMC.h". The header file "BHMC.h" is a class for the tailored basin-hopping Monte Carlo algorithm, where the system parameters are firstly defined, the functions for generating initial configuration of NPs and for calculating the interparticle and interfacial forces are then generated, and the functions associated with conjugate gradient minimization and MC moves are also included. The "main.cpp" performed the operations of this BHMC algorithm. 
2. The code takes an input file "Input", which includes 8 system parameters in total, one per row. As specified below,
* * *
first row: Number of NP I
second row: Number of NP II
third row: Radius of NP I
fourth row: Radius of NP II
fifth row: Omega --- the ratio of the differences in the surface energy of the two types of NPs with the two layers
sixth row: Eff --- factor accounting for the effective radius of the surface-functionalized NPs, which was set to 1.0 throughout the study
seventh row: Chi --- ratio of the difference in the surface energy of NP I with the two fluids to the surface tension between the two fluids
eighth row: Epsilon --- the relative strength of interparticle interactions
* * *
## Code Execution 
* * *
1. Download "main.cpp", "BHMC.h", "Makefile", and "Input" and navigate to the folder including these files
2. In the terminal,
*
$ make
*
which should generate a complied executable code named "BHMC". To run the code, type in the terminal
*
$ ./BHMC Input
*
* * * 
## Output Files
* * * 
The code generates 4 output files:
1. GM_ENERGIES.DAT --- Current best energy minima (till that MC step), where the last one indicates the best energy minima optimized by the algorithm
2. GM_MOVIE.xyz --- NP configurations in xyz format corresponding to the energy minima in the file "GM_ENERGIES.DAT", where the last frame indicates the best NP structure optimized by the algorithm
3. LM_ENERGIES.DAT --- All accepted local energy minima
4. LM_MOVIE.xyz --- NP configurations in xyz format corresponding to the accepted local energy minima in the file "LM_ENERGIES.DAT"
* * *
In addition to the 4 files, the code also prints on the screen the information about initializing NP configuration, conjugate gradient operation and Monte Carlo moves
* * *
## Additional Information
* * *
1. The default box dimension is "15*15*10 times (Radius of NP I)". Periodic boundary conditions were applied to x and y dimensions, and fixed boundary was applied to z dimension. If you need to change the box dimensions (e.g., for larger system size), you can go to Line 53 in "BHMC.h" and change the parameter initialization of "HEIGHT" (for dimension z) and "DENS" (for dimension x and y)
2. An example for searching the hexagonally arranged bilayer is included.
3. Source data for binary NP superlattices discovered in this work is also attached.
