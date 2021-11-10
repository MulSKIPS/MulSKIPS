from dolfin import *
from subprocess import call
import sys, os, time
import numpy as np
import meshio
import cashocs

import pymulskips as pm

#---------------------------------------------------------------
# Provide the name of geo (or directly msh) file
#---------------------------------------------------------------

# structurename = 'Cilinder2_6SiO2_1'
# call([f"gmsh -3 {structurename}.geo -o {structurename}.msh -format msh"], shell=True)

# structurename = 'Si.flat.2micron.KMC240x240x760'
# call([f"gmsh -3 {structurename}.geo -o {structurename}.msh -format msh"], shell=True)

structurename = "TestMundfab"
# structurename = "test_mundfab_armin"



#---------------------------------------------------------------
if structurename == "TestMundfab":
    rescale = 1e6
    mesh, subdomains, boundaries, dx, ds, dS = pm.io.msh2dolfin(structurename+'.msh', save_h5=True, rescale=rescale, rotate_angle=90, rotate_axis=1)
else:
    mesh, subdomains, boundaries, dx, ds, dS = pm.io.msh2dolfin(structurename+'.msh', save_h5=True)



# Setup mulskips process
procdetails = {'substrate': 'Si', 'precursors': ['SiH4', 'H2'], 'temperature': 1098.0536088, 'pressure_SiH4': 6.65, 'pressure_H2': 2666.44}
calibration_params = {
    #Energy barriers (absolute value, in eV if not specified) for MC events:-----------------------------------------
    # Attachment of Si to a site with coor 1, 2 or 3: 
    'Ed1': 3.7-1.74, 'Ed2': 3.7-1.74, 'Ed3': 3.7-1.74, 
    # Detachment of Si from a site without H neighbours with coor 1, 2 or 3:
    'Ee1': 1.4, 'Ee2': 1.74, 'Ee3': 2.4, 
    # Adsorption of H in a site with coor 1, 2 or 3:   
    'Eabs1': 0, 'Eabs2': 0, 'Eabs3': 0,  #NB: these will be neglected for H (Eads=aads-bads*kb*T)!!!
    'aads': 0.8, 
    'bads': 1.5, #non-dimensional
    # Desorption of H from a site with coor 1, 2 or 3:
    'Edes1': 2.48, 'Edes2': 2.48, 'Edes3': 2.48, 
    # Perturbation to the energy barrier (with plus sign) for evaporation of a Si with a single H neighbour: 
    # (for sites with 2 or 3 H neighbours it will be cumulated linearly)
    'deltaH': 0.067, 
    # Dissociation energy barriers (eV) and stoichyometric coefficients (non-dimensional):--------------------------
    'Ereact_SiH4': 0.18, 
    'Ereact_H2': 0, 
    'k_Si_SiH4': 1, #stoichyometric factor k for SiH4->kSi
    'k_H_SiH4': 0, #stoichyometric factor k for SiH4->kH #we're neglecting H contribution from SiH4 (pSiH4<<pH2) to always get the same pads
    'k_H_H2': 1, #stoichyometric factor k for H2->kH
    # Prefactors (non-dimensional if not specified):----------------------------------------------------------------
    'Adep': 1.345068987e9, 
    'Aev': 0.000103*1.0078125, #K^(-1) 
    'Bev': -0.0932*1.0078125,
    's0': 0.16*0.3, 
    'Ades': 2e15*0.3,  
    # Scaling factor (non-dimensional, applied to all events):------------------------------------------------------
    'scalef': 1, 
    }
mpclass = pm.process.CVD(procdetails, calibration_params)



""" Setup dictionary associating each subdomain to a mulskips type 0 (air), 1 (evolving) or 10 (non-evolving/walls)
# Automatically considers only tetras --> subdomains
# Remember that subdomains and boundaries are user-defined. 
# Check them in the info from msh file printed out when running this script (see above) 
"""
subdomains_unique_indices = np.unique(subdomains.array())
print('User-defined subdomains labels in msh file = ', subdomains_unique_indices)
# To make this more general and transparent, we define regions as a dictionary below
regions = { subdomains_unique_indices[0] : 0, # 10, 
            subdomains_unique_indices[1] : 1, 
            subdomains_unique_indices[2] : 1 } # 0 } 
""" 
NB: 0,1,1 makes the pillar disappear immediately (kmc_0-1-1), only a flat surface remains and grows. This is ok for our initial tests.
We'll redo it once the pillar is full 360 deg
"""


# Setup geometry to be used as input in MulSKIPs. Get also the maps, which will be reused in mulskips2dolfin (also to speed it up) 
cell_map, rank_map, wall_map = pm.io.dolfin2mulskips(f'{structurename}_InputKMC.dat', mpclass, regions, mesh, subdomains, regionsfilename=f'{structurename}_KMCregions')


""" Setup & Run """
execpath = os.getenv("PWD")+'/mulskips-source-merge-v2'
rundirname = 'kmc_regions_' + '_'.join([str(x) for x in regions.values()])
endok = pm.setuprun.run_mulskips(execpath, rundirname, \
    Simulation='IN', mp=mpclass, \
    PtransZig=0.9, RunType='R', IDUM=9117116, \
    ExitStrategy='Iter', OutMolMol=2500000, IterMax=5000000, \
    cadfilename=f'{structurename}_InputKMC.dat', \
    SaveCoo=True, coofilename=f'{structurename}_AfterKMC.dat', \
    SaveFinalState=False)
if not endok:
    print('\nSTOP!!! Few MC particles in KMC...')
    sys.exit()


""" Update mesh with new geometry obtained with MulSKIPs
- mulskips2dolfin() should be exactly the same as read_phases_from_mulskips()
  In this way we will update a meshfunction (e.g., Ph_tmp in LA) leaving the mesh unaltered (walls will be mapped from LattFunc_array)
- Leave the mesh unaltered and write the msh with the same meshcoordinates. Synopsys will do a remeshing by themselves. 
  This could help keeping coherence in the syn workflow...talk to Eberhard
- Otherwise a solution could be to try to refine ourselves the mesh around the new boundaries, based adding new "solid" points around the current new solid coarse meshpoints
      and refining the mesh nearby them...Check this out:
      https://fenicsproject.discourse.group/t/why-are-boundary-and-surface-markers-not-carried-over-to-the-refined-mesh/5822/2
"""
subdomains = pm.io.mulskips2dolfin(f'{structurename}_AfterKMC.dat', mpclass, mesh, subdomains, regions, \
    cell_map, rank_map, wall_map, pvdfilename=f'{structurename}_Mesh_AfterKMC', LA=False)

""" Write dolfin mesh to msh format. 
This updates the old mesh with the new subdomains.
Writes the updated subdomains to pvd (look at it in paraview).
Writes the updated mesh to msh (only subdomains markers have changed!)
NB: remember to rotate/scale it again at the end if needed
"""
pm.io.dolfin2msh(structurename+'.msh', structurename+'_final.msh', \
    mesh, subdomains, boundaries)

