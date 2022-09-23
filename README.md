### MulSKIPS
A Kinetic Monte Carlo super-Lattice code, designed to simulate with atomic resolution the kinetics of processes (e.g., PVD, CVD, laser annealing) involving elements, alloys and compounds characterized by the sp3 bond symmetry

### INSTALLATION 

- To compile MulSKIPS use the following commands:

```
git clone https://github.com/MulSKIPS/MulSKIPS.git
cd /your_MulSKIPS_directory/mulskips-source/ 
make clean ; make 
```

If needed, edit the makefile.
A successful compilation will generate a "mulskips.e" executable file.

- To run a simulation, go to your working directory, use pymulskips routines to generate a "start.dat" file with all simulation parameters, then simply call  /your_MulSKIPS_directory/mulskips.e  

- To use the pymulskips routines in your python script, simply update your PYTHONPATH:

```
export $PYTHONPATH=$PYTHONPATH:/your_MulSKIPS_directory/pymulskips
```

- For examples of usage, see /your_MulSKIPS_directory/examples/.




### MORE DOCUMENTATION and EXAMPLE on PVD, CVD, LASER ANNEALING processes for various materials (Si, SiC, SiGe, ...) will be uploaded soon, so...STAY TUNED!


### REFERENCES

[A. La Magna, A. Alberti, E. Barbagiovanni, C. Bongiorno, M. Cascio, 
I. Deretzis, F. La Via, and E. Smecca, "Simulation of the growth kinetics
in group iv compound semiconductors" Phys. Status Solidi A 216, 1800597 (2019)].

[G. Calogero, D. Raciti, P. Acosta-Alba, F. Cristiano, I. Deretzis, G. Fisicaro, K. Huet,
S. Kerdilès, A. Sciuto and A. La Magna, "Multiscale modeling of ultrafast melting
phenomena", npj Computational Materials 8, 36 (2022)]


<!--
block comment
-->
