### MulSKIPS
A Kinetic Monte Carlo super-Lattice code, designed to simulate with atomic resolution the kinetics of processes (e.g., PVD, CVD, laser annealing) involving elements, alloys and compounds characterized by the sp3 bond symmetry

INSTALLATION 

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


<!--
**MulSKIPS/MulSKIPS** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
