### MulSKIPS
A Kinetic Monte Carlo super-Lattice code, designed to simulate with atomic resolution the kinetics of processes (e.g., PVD, CVD, laser annealing) involving elements, alloys and compounds characterized by the sp3 bond symmetry

INSTALLATION

To compile MulSKIPS use the following commands:
- git clone <URL> (clone this directory)
- cd /your_MulSKIPS_directory/mulskips-source/ (go to "mulskips-source" directory)
- (edit makefile if needed)
- make clean ; make 
- cd ..
This will create the "mulskips.e" executable file.

To use the routines in pymulskips simply add its path to your PYTHONPATH:
- export $PYTHONPATH=$PYTHONPATH:/your_MulSKIPS_directory/pymulskips
or copy the directory within the working directory from where you are running mulskips.e

For example of usage, go to /your_MulSKIPS_directory/examples/

MORE DOCUMENTATION and EXAMPLE on PVD, CVD, LASER ANNEALING processes for various materials (Si, SiC, SiGe, ...) will be uploaded soon, so...STAY TUNED!


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
