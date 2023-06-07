from dolfin import *
from subprocess import call
import sys, os, time
import numpy as np
from pymulskips import io, setuprun, process, analyze
import cantera
print(f'cantera v{cantera.__version__}')
import meshio
print(f'meshio v{meshio.__version__}')
import cashocs
print(f'cashocs v{cashocs.__version__}')


# Set the path of MulSKIPS source code, and test if it exists.
execpath = './mulskips-source'
if not os.path.exists(execpath):
    print(f"Error in run.py, the path of mulskips-source folder: {execpath} does not exist.")
    sys.exit(1)


""" 
Setup & Run MulSKIPS 
Quickly Recompiles MulSKIPS to reset the box size.  
Then runs the MulSKIPS simulation.
All results will be in rundirname.
"""

# NB: with mp='custom' we are telling MulSKIPS to read the PRE-EXISTING start.dat placed in rundirname folder!

setuprun.setup_mulskips_src(execpath, lenx=780, leny=780, lenz=780) # it will rerecompile if needed
rundirname = '.'
Nout = 50
itermax = 1200000000 
endok = setuprun.run_mulskips(execpath, rundirname, \
    Simulation='SG', mp='custom', Seed_box=[48,0,0], \
    PtransZig=0.95, RunType='R', IDUM=9117116, \
    ExitStrategy='Iter', OutMolMol=int(itermax/Nout), IterMax=itermax, \
    SaveFinalState=False, setup_only=False)
if not endok:
    print('\nSTOP!!! Few MC particles in KMC...')
    sys.exit()

