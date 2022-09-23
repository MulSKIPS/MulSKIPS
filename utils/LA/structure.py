from dolfin import *
from subprocess import call
import math
import numpy
import cmath

lambda0=1.0
epsilon_sphere=complex(-10.0, -1.0)
n_bkg=1
kdir_inc=(-1, 0, 0) #propagation direction incident field
Edir_inc=Constant((0.0, 0.0, -1.0)) #polarization incident field

k0=2*pi/lambda0
k2=n_bkg*k0



# Load mesh and define function space
call(["gmsh -3 flat.large.geo -o structure.msh -format msh2"],shell=True)
call(["dolfin-convert structure.msh structure.xml"],shell=True)

mesh = Mesh("structure.xml")
mv=MeshFunction("size_t",mesh,"structure_physical_region.xml")
#ms=MeshFunction("size_t",mesh,"structure_facet_region.xml")
hdf = HDF5File(mesh.mpi_comm(), "file.h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(mv, "/mv")
#hdf.write(ms, "/ms")
hdf.close()

