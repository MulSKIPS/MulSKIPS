This folder contains different examples of laser annealing simulations performed with MulSKIPS.

In particular, each folder contains the minimal code necessary to reproduce the results reported in this recent article (submitted):
"Multiscale modelling of ultrafast melting and structural disorder in group IV alloys" 
by Gaetano Calogero, Domenica Raciti, Damiano Ricciarelli, Pablo Acosta-Alba,
Fuccio Cristiano, Richard Daubriac, Remi Demoulin, Ioannis Deretzis, Giuseppe Fisicaro,
Jean-Michel Hartmann, Sébastien Kerdilés, Antonino La Magna 

### USAGE

To run the simulations, you should:
- Make sure you installed all the dependencies (see below)
- Make sure you copied the pymulskips directory in this folder, or added it to your PYTHONPATH
- Make sure you updated the "execpath" variable in run.py with your mulskips-source directory FULL path

- Use the following command to construct the msh file from the geo file located in each folder:
```
python3 structure.py
```

- Use the following command to run the simulation as background process (note that you need to ensure a few GBs of free storage to run this and the other examples):
```
nohup python3 -u MulSKIPS_LA_sockets.py &
```
- The script requires hours, so go grap a coffee :)

- At the end of the simulations you will have a sequence of KMC_t_* folders with XYZ files inside. Coordinates of cubic undercoordinated atoms are inside I00000000.xyz like files, whereas hexagonal bulk sites (only if PtransZig was set < 1 in the python script!) are in the I0000000_d.xyz files.
- PVD and VTU files contain mesh functions (e.g., InterpXGeS* files contain the local Ge fraction in the FEM mesh).
- Use the following command to run the analysis script, which will extract relevant physical info from log.out and log.mulskips.out text files:
```
./plot.sh
```
 
- The code does not support MPI.



### DEPENDENCIES

This example was tested for Python 3.8.10 using Ubuntu.

The libraries in pymulskips/requirements_aptdolfin.txt should be used with dolfin 2019.2.0.dev0 version.

The libraries in pymulskips/requirements_condadolfin.txt should be used with dolfin installed through CONDA.

The versions indicated are the only ones tested. 

Please feel free to try some newest ones and report an issue if something does not work.

- To install FEniCS on Unix system (tested also on WSL2 for Windows 10/11, with Ubuntu 20.04):

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
```

Ref. https://fenicsproject.org/download/archive/


- To install Cantera:

```
pip3 install cantera
```

Ref. https://cantera.org/install/


- To install cashocs:

```
pip3 install cashocs
```

Make sure to install all cashocs dependencies first
(see https://pypi.org/project/cashocs/#manual-installation)

