from pymulskips_sockets import analyze
import os

rundir = os.getenv('PWD')+'/kmc_regions_10_1_10_0_10_10_10_10'
Niter = 10
xyzfile = rundir+f'/I000000{Niter}.xyz'
alat = 5.43 # Ang
what = 'surface+coverage'
# what = 'coverage'
# what = 'surface'
DEP3Dfile = rundir+f'/I000000{Niter}.{what}.DEP3D.xyz'
newfile = rundir+f'/I000000{Niter}.{what}.xyz'

analyze.export_xyz(xyzfile, newfile, alat, what=what, DEP3Dfile=DEP3Dfile)
