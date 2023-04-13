## [MulSKIPS](https://mulskips.readthedocs.io/en/latest/index.html) ##
A Kinetic Monte Carlo super-Lattice code, designed to simulate with atomic resolution the kinetics of processes (e.g., PVD, CVD, laser annealing) involving elements, alloys and compounds characterized by the sp3 bond symmetry

## Installation ##

- To compile MulSKIPS use the following commands:

```
git clone https://github.com/MulSKIPS/MulSKIPS.git
cd /your_MulSKIPS_directory/mulskips-source/ 
make 
```

If something goes wrong, please double check (and edit, if needed) the makefile, and ensure that your GNU Fortran compiler (gfortran or f95) is updated (tested version >= 9.4.0).
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
- For examples of usage, see /your_MulSKIPS_directory/examples/.


<!--
## Tutorials and examples ##
The easiest way to get started is to follow the tutorials [here](https://mulskips.readthedocs.io/en/latest/index.html).
-->

## Documentation ##

Please find documentation here:

- [Documentation](https://mulskips.readthedocs.io/en/latest/index.html)


## References ##

[A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, 
I. Deretzis, F. La Via, and E. Smecca, "Simulation of the growth kinetics
in group iv compound semiconductors" Phys. Status Solidi A 216, 1800597 (2019)].

[G. Calogero, D. Raciti, P. Acosta-Alba, F. Cristiano, I. Deretzis, G. Fisicaro, K. Huet,
S. Kerdilès, A. Sciuto and A. La Magna, "Multiscale modeling of ultrafast melting
phenomena", npj Computational Materials 8, 36 (2022)]

[G. Calogero, D. Raciti, D. Ricciarelli, P. Acosta-Alba, F. Cristiano, R. Daubriac, 
R. Demoulin, I. Deretzis, G. Fisicaro, J.-M. Hartmann, S. Kerdilés, A. La Magna, 
"Multiscale modelling of ultrafast melting and structural disorder in group IV alloys" (submitted)]

<!--
block comment
-->
