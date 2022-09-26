### USAGE

- Make sure you installed all the dependencies below
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

- To run this example you need at least 5 GB free storage, to allocate 2 large auxiliary DAT files which can be removed after the run if you want.
These contain the geometry information needed to import the mesh into MulSKIPS.
BUT keep them in the folder to avoid repeating the slow (a couple of hours on a good machine) MSH-MulSKIPS interpolation 
in case you want to repeat only the KMC run, which is mush faster (minutes).
Please note that the code does not support MPI.


### DEPENDENCIES

Tested only for Python 3.8.10 using Ubuntu.
The following libraries are required:
- dolfin (2019.2.0.dev0)
- cantera (2.6.0a4)
- cashocs (>= 1.4.0)

The versions indicated are the only ones tested. 
Please feel free to try some newest ones and report an issue if something does not work.

NB: Conda installation of libraries WAS NOT TESTED.

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


