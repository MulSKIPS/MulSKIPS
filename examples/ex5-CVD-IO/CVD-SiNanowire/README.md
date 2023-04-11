### USAGE

- Make sure you installed all the dependencies (see below)
- Make sure you copied the pymulskips directory in this folder, or added it to your PYTHONPATH
- Make sure you updated the "execpath" variable in run.py with your mulskips-source directory FULL path
- Use the following command:

```
python3 run.py
```

or this to run in as background process

```
nohup python3 -u run.py > log.out 2>&1 &
```

- The script produces a sequence of XYZ files located within the output /kmc_regions_\*/ directory (the number of snapshots "Nout" and can be changed withi the run.py script). 
- A \*DEP3D.xyz file is also produced within the /kmc_regions_\*/ directory, ready to be imported in the DEP3D tool.  
- To run this example you need at least 5 GB free storage, to allocate 2 large auxiliary DAT files which can be removed after the run if you want. These contain the geometry information needed to import the mesh into MulSKIPS.
BUT keep them in the folder to avoid repeating the slow (a couple of hours on a good machine) MSH-MulSKIPS interpolation in case you want to repeat only the KMC run, which is mush faster (minutes).
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


