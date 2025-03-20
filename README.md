## [MulSKIPS](https://mulskips.readthedocs.io/en/latest/index.html) ##
A Kinetic Monte Carlo super-Lattice code, designed to simulate with atomic resolution the kinetics of processes (e.g., PVD, CVD, laser annealing) involving elements, alloys and compounds characterized by the sp3 bond symmetry. It can simulate generation and evolution of point and extended defects (e.g. stacking faults), as well as the symultaneous evolution of multiple crystalline phases (e.g. cubic and hexagonal) during the process. Basic geometries (surfaces, nanocrystals) or complex TCAD meshes can be used as input. Setup and post-processing is managed by the user-friendly Python module ```pymulskips``` (see examples in the repository).

## Installation ##

- To compile MulSKIPS use the following commands:

```
git clone https://github.com/MulSKIPS/MulSKIPS.git
cd /your_MulSKIPS_directory/mulskips-source/ 
make 
```

If something goes wrong, please double check (and edit, if needed) the makefile, and ensure that your GNU Fortran compiler (gfortran or f95) is updated (tested version >= 9.4.0). Please note that if you are using an arm64 architecture you will need to replace -mcmodel=medium with -mcmodel=small in the FFLAGS line of the makefile.  
A successful compilation will generate a "mulskips.e" executable file in the /your_MulSKIPS_directory/mulskips-source/ directory.

- To run a simulation, go to your working directory, use pymulskips routines to generate a "start.dat" file with all simulation parameters, then simply call  /your_MulSKIPS_directory/mulskips.e  
- To use the pymulskips routines in your python script, simply copy the directory /your_MulSKIPS_directory/pymulskips/ into your working directory.
In alternative, you can add /your_MulSKIPS_directory/pymulskips/ to your PYTHONPATH environment variable.
- The pymulskips module needs a lot of dependencies. 
The required versions of these dependencies may vary with the installed version of dolfin.
If you installed dolfin through CONDA, please note that the following dependencies were successfully tested (install them using the following command):
```
pip3 install -r /your_MulSKIPS_directory/pymulskips/requirements_condadolfin.txt
```
If you installed dolfin through apt, please note that dolfin 2019.2.0.dev0 version has been successfully tested using the following versions for the dependencies (install them using the following command):
```
pip3 install -r /your_MulSKIPS_directory/pymulskips/requirements_aptdolfin.txt
```

## Tutorials and examples ##

Check [here](https://github.com/MulSKIPS/MulSKIPS/tree/main/examples) for some examples of usage. 

<!--
## Tutorials and examples ##
The easiest way to get started is to follow the tutorials [here](https://mulskips.readthedocs.io/en/latest/index.html).
-->

## Documentation ##

Please find documentation [here](https://mulskips.readthedocs.io/en/latest/index.html).


## References ##

[1] A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, I. Deretzis, F. La Via, and E. Smecca, "Simulation of the Growth Kinetics in Group IV Compound Semiconductors", physica status solidi (a) vol. 216, no. 10, p. 1800597, 2019, doi: 10.1002/pssa.201800597

[2] G. Fisicaro, C. Bongiorno, I. Deretzis, F. Giannazzo, F. La Via, F. Roccaforte, M. Zielinski, M. Zimbone, A. La Magna, "Genesis and Evolution of Extended Defects: The Role of Evolving Interface Instabilities in Cubic SiC", Applied Physics Reviews vol. 7, no. 2, p. 021402, Apr. 2020, doi: 10.1063/1.5132300

[3] https://hq.imm.cnr.it/content/super-lattice-kinetic-monte-carlo-method-simulate-cvd-epitaxy-si-based-materials

[4] G. Calogero, D. Raciti, P. Acosta-Alba, F. Cristiano, I. Deretzis, G. Fisicaro, K. Huet, S. Kerdilès, A. Sciuto and A. La Magna, "Multiscale modeling of ultrafast melting phenomena", npj Computational Materials 8, 36 (2022), doi: 10.1038/s41524-022-00720-y

[5] G. Calogero, D. Raciti, D. Ricciarelli, P. Acosta-Alba, F. Cristiano, R. Daubriac, R. Demoulin, I. Deretzis, G. Fisicaro, J.-M. Hartmann, S. Kerdilés, A. La Magna, "Atomistic insights into ultrafast SiGe nanoprocessing" (submitted)

<!--
block comment
-->
