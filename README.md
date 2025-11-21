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
Note that dolfin is required only to use the io.py module for interactions with meshes and PDE solvers (e.g., for laser annealing). 
In fact, dolfin will be searched for only if you try to import the pymulskips.io module.
If you choose to install dolfin, you should install the legacy version 2019.1.0, which was released in April 2019:
https://fenicsproject.org/download/archive/
This is because some of the MulSKIPS' dependencies do not work with dolfinx.
The required versions of the other dependencies may vary with the installed version of dolfin.
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

## Developers ##
Main developers:
- Gaetano Calogero (CNR-IMM)
- Antonino La Magna (CNR-IMM)
- Giuseppe Fisicaro (CNR-IMM)

Contributors:
- Damiano Ricciarelli (CNR-IMM)
- Ioannis Deretzis (CNR-IMM)
- Domenica Raciti
- Remi Helleboid

CNR-IMM: Istituto per la Microelettronica e Microsistemi, Consiglio Nazionale delle Ricerche, Catania, Italy.

## Publications citing ```MulSKIPS``` ##

- A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, I. Deretzis, F. La Via, and E. Smecca, "Simulation of the Growth Kinetics in Group IV Compound Semiconductors", physica status solidi (a) vol. 216, no. 10, p. 1800597, 2019, DOI: https://doi.org/10.1002/pssa.201800597

- G. Fisicaro, C. Bongiorno, I. Deretzis, F. Giannazzo, F. La Via, F. Roccaforte, M. Zielinski, M. Zimbone, A. La Magna, "Genesis and Evolution of Extended Defects: The Role of Evolving Interface Instabilities in Cubic SiC", Applied Physics Reviews vol. 7, no. 2, p. 021402, Apr. 2020, DOI: https://doi.org/10.1063/1.5132300

- G. Calogero, D. Raciti, P. Acosta-Alba, F. Cristiano, I. Deretzis, G. Fisicaro, K. Huet, S. Kerdilès, A. Sciuto and A. La Magna, "Multiscale modeling of ultrafast melting phenomena", npj Computational Materials 8, 36 (2022), doi: https://doi.org/10.1038/s41524-022-00720-y

- G. Calogero, I. Deretzis, G. Fisicaro, M. Kollmuß, F. La Via, S.F. Lombardo, M. Sch\"oler, P. Wellmann, A. La Magna, "Multiscale simulations for defect-controlled processing of group IV materials", Crystals 12 (12), 1701 (2022) DOI: https://doi.org/10.3390/cryst12121701 

- G. Calogero, D. Raciti, D. Ricciarelli, P. Acosta-Alba, F. Cristiano, R. Daubriac, R. Demoulin, I. Deretzis, G. Fisicaro, J.-M. Hartmann, S. Kerdil`es, and A. La Magna, "Atomistic insights into ultrafast SiGe nanoprocessing", The Journal of Physical Chemistry C, 127 (39), 19867 (2023) DOI: https://doi.org/10.1021/acs.jpcc.3c05999 

- D. Raciti, G. Calogero, D. Ricciarelli, R. Anzalone, G. Morale, D. Murabito, I. Deretzis, G. Fisicaro, A. La Magna, "Multiscale atomistic modelling of CVD: from gas-phase reactions to lattice defects", Materials Science in Semiconductor Processing, 167, 107792 (2023) DOI: https://doi.org/10.1016/j.mssp.2023.107792 

- G. Calogero, I. Deretzis, G. Fisicaro, D. Ricciarelli, R.G. Viglione, A. La Magna, "Tailoring nuclear spins order with defects: a Quantum Technology CAD study", Advanced Quantum Technologies, e2500160 (2025) DOI: https://doi.org/10.1002/qute.202500160

<!--
block comment
-->
