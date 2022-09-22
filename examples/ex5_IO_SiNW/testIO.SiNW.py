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


""" 
Load mesh and import it in Dolfin
"""

structurename = 'si-cylinder-ox-removed-withphysicalentitiers_fps' # neglect .msh extension

# Use the line below to convert a geo into a msh file
# call([f"gmsh -3 structurename.geo -o structurename.msh -format msh"], shell=True)

rescale = 1e3 # msh by Chiara is in micron. Here we need nm. 
mesh, subdomains, boundaries, dx, ds, dS = io.msh2dolfin(structurename+'.msh', \
    save_h5=True, rescale=rescale, rotate_angle=90, rotate_axis=1)



""" 
Setup MulSKIPS process
""" 

process_data = {'substrate': 'Si', \
                'precursors': ['H2', 'HCl', 'SiCl2'], \
                'temperature': 1253.15, \
                'partial_pressures': {'H2': 7813.942439181746, 'HCl': 154.92481855016135, 'SiCl2': 28.731950245466393}}
calibration = {
# Energy barrier (absolute value) for deposition of Si in a site with coor 1, 2 or 3:
'Edep': {'Si': [0]*3}, # eV
# Energy barrier (absolute value) for evaporation of Si from a site without H neighbours with coor 1, 2 or 3:
'Eev': {'Si': np.asarray([-0.4,  0.3,  1.1])}, # eV
# Energy barrier (absolute value) for absorption of H in a site with coor 1, 2 or 3:
'Eabs': {'Cl': [0]*3, 'H': [0]*3}, # eV  
#Energy barrier (with minus sign) for desorption of H from a site with coor 1, 2 or 3:
'Edes': {'Cl': np.asarray([0]*3), 'H': np.asarray([0]*3)}, # eV
# Perturbation to the energy barrier for evaporation of a Si with coverage neighbours
'deltaE': {'H': np.asarray([0]*3), 'Cl': np.asarray([0]*3)}, # eV
# Parameters for s0, used in attachment and absorption (dict, one entry per molecule)
'alpha': {'SiCl2': 2.6193843228288402e+23, 'H2': 9.666800569359601e+21, 'HCl': 2.71531702513564e+22}, #non-dimensional
'kd0kr0ratio': {'SiCl2': 600, 'H2': 1200000000, 'HCl': 8400000}, #non-dimensional
'EdminusEr': {'SiCl2': 0.1647809, 'H2': 1.57, 'HCl': 1.6}, #eV
# Parameters for kevap, used in detachment and desorption (dict, one entry per molecule)
'kevap_A': {'H2': 2.4767826e+32, 'HCl': 2.3486188e+34, 'SiCl2': 1.0360509e+32}, #[m^2/s]  
'kevap_E': {'H2': 2.69, 'HCl': 3.12, 'SiCl2': 2.91}, #[eV]
# Global scaling factor, used to scale all events probabilities
'scalef': 1.597677326882467e-19, # non-dimensional 
'scalefcov': 1.597677326882467e-19   # non-dimensional 
}
mygas= cantera.Solution('CVD_ST.xml')
mpclass = process.CVD(process_data, calibration_params=calibration, gas=mygas)


""" 
Setup dictionary associating each subdomain to a MulSKIPS region type 
>>> 0 (gas)
>>> 1 (evolving regions) 
>>> 10 (non-evolving regions)
This automatically considers only tetrahedrons in msh (subdomains).
Visualize them (e.g. in gmsh) to get the indices of the regions in the msh to be used below.
"""

subdomains_unique_indices = np.unique(subdomains.array())
print('User-defined subdomains labels in msh file = ', subdomains_unique_indices)
# To make this more general and transparent, we define regions as a dictionary below
# Remember to assign a type to all tetra regions in the MSH (visualize them in SVisual or gmsh)
regions = { subdomains_unique_indices[0] : 0, # gas
            subdomains_unique_indices[1] : 1} # evolving silicon

""" 
Convert structure from Dolfin to MulSKIPS readable structure (outputs PVD/VTU files, and a DAT file).
Gets also some auxiliary maps, which can be reused in mulskips2dolfin if needed (not needed when using DEP3D)
Returns also the box size to be used in MulSKIPS. 
NB: since this takes time, recycle existing DAT file in the current directory
"""

if not os.path.isfile(os.getenv("PWD")+f'/{structurename}_InputKMC.dat'):
    cell_map, rank_map, wall_map, lenx, leny, lenz = io.dolfin2mulskips(f'{structurename}_InputKMC.dat', \
        mpclass, regions, mesh, subdomains, \
        regionsfilename=f'{structurename}_KMCregions', return_box_size=True)
else:
    print("Found DAT input file for MulSKIPS in current folder...")
    alat = mpclass.KMC_sf * 1e-1 # nm
    lenx, leny, lenz = io.get_boxsize_from_mesh(mesh, alat)


""" 
Setup & Run MulSKIPS 
Quickly Recompiles MulSKIPS to reset the box size.  
Then runs the MulSKIPS simulation.
All results will be in rundirname.
"""

execpath = os.getenv("PWD")+'/mulskips-source' # path of MulSKIPS source code 
setuprun.setup_mulskips_src(execpath, lenx, leny, lenz) # it will rerecompile if needed
rundirname = 'kmc_regions_' + '_'.join([str(x) for x in regions.values()])
Nout = 10 # number of outputs
itermax = 350000000 # total number of KMC iterations
endok = setuprun.run_mulskips(execpath, rundirname, \
    Simulation='IN', mp=mpclass, \
    PtransZig=1.0, RunType='R', IDUM=9117116, \
    ExitStrategy='Iter', OutMolMol=int(itermax/Nout), IterMax=itermax, \
    cadfilename=f'{structurename}_InputKMC.dat', \
    SaveCoo=True, coofilename=f'{structurename}_AfterKMC.dat', \
    SaveFinalState=False, setup_only=False)
if not endok:
    print('\nSTOP!!! Few MC particles in KMC...')
    sys.exit()



"""
Extract superficial atoms only from the final structure, to be used in DEP3D
"""

Niter = Nout
rundir = os.getenv('PWD')+'/'+rundirname
xyzfile = rundir+f'/I000000{Niter}.xyz'
what = 'surface+coverage'
# what = 'coverage'
# what = 'surface'
DEP3Dfile = rundir+f'/I000000{Niter}.{what}.DEP3D.xyz'
newfile = rundir+f'/I000000{Niter}.{what}.xyz'
analyze.export_xyz(xyzfile, newfile, alat=mpclass.KMC_sf, what=what, DEP3Dfile=DEP3Dfile)





""" OPTIONAL 
NB: dolfin2mulskips() MUST be run before this 
Update mesh with new geometry obtained with MulSKIPs
Leave the mesh unaltered and write the final MSH with the same mesh coordinates. 
Synopsys will then need to do a remesh. 
"""

subdomains = io.mulskips2dolfin(f'{structurename}_AfterKMC.dat', mpclass, mesh, subdomains, regions, \
    cell_map, rank_map, wall_map, pvdfilename=f'{structurename}_Mesh_AfterKMC', LA=False)

""" OPTIONAL
Write dolfin mesh to msh format. 
Makes new mesh with updated physical entities
"""

io.dolfin2msh(structurename+'.msh', structurename+'_final.msh', \
    mesh, subdomains, boundaries)



