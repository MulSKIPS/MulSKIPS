### USAGE

- Make sure you installed all the dependencies below
- Make sure you updated the "execpath" in testIO.FDSOI.py with your mulskips-source directory
- Use the following command:

```
python3 testIO.FDSOI.py
```

OR

```
nohup python3 -u testIO.FDSOI.py > log.out 2>&1 &
```

- To run this example you need at least 5 GB free storage, to allocate 2 large auxiliary DAT files which can be removed after the run if you want.
These contain the geometry information needed to import the mesh into MulSKIPS.
BUT keep them in the folder to avoid repeating the slow (a couple of hours on a good machine) MSH-MulSKIPS interpolation 
in case you want to repeat only the KMC run, which is mush faster (minutes).
No MPI is used.


### DEPENDENCIES

Tested only for Python 3.8.10 using Ubuntu with the following libraries:
- dolfin (2019.2.0.dev0)
- cantera (2.6.0a4)
- cashocs (1.4.0)
- meshio (5.0.2)

The versions indicated are the only ones tested. 
Please feel free to try some newest ones and report an issue if something does not work.

NB: Conda installation of libraries WAS NOT TESTED.

To install FEniCS on Unix system (tested also on WSL2 for Windows, with Ubuntu 20.04):
(ref. https://fenicsproject.org/download/archive/)
To install FEniCS on Ubuntu, run the following commands:

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
```

For meshio and cashocs:
(ref. https://pypi.org/project/cashocs/#manual-installation)

```
pip3 install meshio[all] --no-binary=h5py
pip3 install cashocs
```

For Cantera:
(ref. https://cantera.org/install/)

```
pip3 cantera
```

