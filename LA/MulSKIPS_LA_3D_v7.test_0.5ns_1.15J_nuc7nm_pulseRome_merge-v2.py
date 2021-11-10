import sys
import os
import shutil
import xml.etree.ElementTree
from dolfin import *
import time
import subprocess
from scipy import linalg
import math
import numpy
import cmath
from ufl import tanh
import time
# from pymulskips.mulskips_pylib import run_mulskips, Si_melting_MulskipsProcess
from pymulskips.mulskips_pylib_merge_v3 import run_mulskips, Si_melting_MulskipsProcess
from pymulskips.LaserAnnealing_fast import *

t000 = time.time()

comm = MPI.comm_world
rank = MPI.rank(comm)
size_mpi= MPI.size(comm)

# Constants: global variables
scale = 1.e9
geometric_tol = 1e-3
#tic = 0.
#if rank == 0 :
#   tic = time.clock()

# Physics Constant (not material dependent)
kBcost = 8.617385e-5  # [eV K^-1]
c_light= 299792458   # m/s
Ph_ini = 1.
a1 = 0.8839          # Karma Constant 1
a2 = 0.3982          # Karma Constant 2
eps0 =8.8541878128e-12 # F m-1
eC = 1.602176634e-19 # electron charge
T_ini = 300.
T_back = 300.
T_resistor = True
AddedSubstrate =1.e5

wave_lenght = 308.
k0_wnumb = 2*pi/wave_lenght
k0sq = k0_wnumb*k0_wnumb
dt_maxT = 1.0 # nanoseconds
dt_min = 0.5 # nanoseconds
t = 0.
t_stop = 10.e2
toffset = 0.
end_simulation = False
verb_log = True
Elas = 11500

epsSiR_c = "-9.025e-7 * pow(T,2) - 6.55813e-3 * T + 13.8925" # 11.942
epsSiI_c = "9.652E-3 * T + 35.0688" # 37.435
epsSiR_l = "-15.734"
epsSiI_l = "10.126"
rh0Si_c = "2320"  # kg/m^3
rh0Si_l = "2520"  # kg/m^3
kSi_c   = "100 * ( ( 1523.7 * pow(T,-1.226) ) * ( 0.5*(1.0+tanh((1200-T)/slp))) + ( 9 * pow(T,-0.502)) * ( 0.5*(1.0+tanh((T-1200)/slp))))" # Average 880  J/kg K
kSi_l = "100 * ( 0.502 + 0.000293 * ( T - Tm ))" # W /mK
cpSi_c  = "10 * pow(T,1.034) / ( 1.02 + 0.01 * T ) - 213" # Average 82 , W /mK
cpSi_l = "1045" # J/kg K
LHeatSi_c = "1797000" # J/kg

# # Read data from database (these expressions are written in C)
# database = xml.etree.ElementTree.parse('Materials.xml').getroot()
# ns = {'material': 'http://www.screen-lasse.com/LIAB/Material'}
# for material in database.findall('material:Material', ns):
#     if material.attrib['Name'] == 'Si':
#         # Real permittivity crystal
#         for prop in material.find('material:Crystal', ns).findall('material:RealPermittivity', ns):
#             if prop.attrib['Name'] == '308':
#                 epsSiR_c = prop.find('material:Formula', ns).attrib['Value']
#                 break
#         # Imaginary permittivity crystal
#         for prop in material.find('material:Crystal', ns).findall('material:ImaginaryPermittivity', ns):
#             if prop.attrib['Name'] == '308':
#                 epsSiI_c = prop.find('material:Formula', ns).attrib['Value']
#                 break

#         # Real permittivity liquid
#         for prop in material.find('material:Liquid', ns).findall('material:RealPermittivity', ns):
#             if prop.attrib['Name'] == '308':
#                 epsSiR_l = prop.find('material:Formula', ns).attrib['Value']
#                 break
#         # Imaginary permittivity liquid
#         for prop in material.find('material:Liquid', ns).findall('material:ImaginaryPermittivity', ns):
#             if prop.attrib['Name'] == '308':
#                 epsSiI_l = prop.find('material:Formula', ns).attrib['Value']
#                 break

linearsolver = 'Direct'
abstol = 1E-5
reltol = 1E-3
divlim = 1E4
no0ini = False

# Read data for Laser Pulse (tabulated pulse)
pulse_name = 'pulse_rome.xml'
# pulse_name = 'ARC160ns_NEW.xml'
tree = xml.etree.ElementTree.parse(pulse_name)
root = tree.getroot()
p_tab=[]
for time_int in root:
    p_tab.append((scale*float(time_int.attrib['Coordinate']),float(time_int.attrib['Value'])))
# rescale time scale pulse


parameters["ghost_mode"] = "shared_vertex"
parameters["form_compiler"]["quadrature_degree"] = 5
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "file.h5", "r")
hdf.read(mesh, "/mesh", False)
mv = MeshFunction("size_t", mesh, mesh.topology().dim()) #CellFunction
hdf.read(mv, "/mv")
hdf.close()
num_lay=2 # 0 is the Air
degree_L=1
dimen = 3
lagr_el = FiniteElement('CG', mesh.ufl_cell(), degree_L)
V = FunctionSpace(mesh, lagr_el) #c
dofmap = V.dofmap()
my_first, my_last = dofmap.ownership_range() # global
coord_mesh = V.tabulate_dof_coordinates().reshape((-1, dimen))
unowned = dofmap.local_to_global_unowned()
dofs = list(filter(lambda dof: dofmap.local_to_global_index(dof) not in unowned, range(my_last-my_first)))
coord_mesh=coord_mesh[dofs]
n_dof = len(coord_mesh)
print("number degree of freedom", n_dof)
print("number cells", mesh.num_cells())

dx = Measure("dx", domain=mesh, subdomain_data=mv)

#dS = Measure('dS')(subdomain_data=boundaries)
# Create classes for defining parts of the boundaries and the interior of the domain
xBox0=0
xBox1=21.72
yBox0=0
yBox1=21.72
if (dimen > 2) :
   zBox0=-210.
   zBox1=19790.390
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], xBox0)
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], xBox1)
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[1] < (yBox0 + 1.e-3) and x[1] > (yBox0 - 1.e-3) and on_boundary)
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[1] < (yBox1 + 1.e-3) and x[1] > (yBox1 -1.e-3) and on_boundary)
if (dimen > 2):
  class Zmin(SubDomain):
      def inside(self, x, on_boundary):
          return bool(x[2] < (zBox0 + 1.e-3) and x[2] > (zBox0 -1.e-3) and on_boundary)
  class Zmax(SubDomain):
      def inside(self, x, on_boundary):
          return bool(x[2] < (zBox1 + 1.e-3) and x[2] > (zBox1 -1.e-3) and on_boundary)

# Initialize boundary sub-domain instances
left = Left()
right = Right()
top = Top()
bottom = Bottom()
if (dimen > 2):
   zmin_bnd = Zmin()
   zmax_bnd = Zmax()
# Initialize mesh function for boundary domains
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1) # FacetFunction("size_t", mesh)
ind_bottom =1
ind_top = 2
ind_left = 3
ind_right = 4
ind_zmin = 5
ind_zmax = 6

boundaries.set_all(0)
left.mark(boundaries, ind_left)
right.mark(boundaries, ind_right)
top.mark(boundaries, ind_top)
bottom.mark(boundaries, ind_bottom)
if (dimen > 2):
   zmin_bnd.mark(boundaries, ind_zmin)
   zmax_bnd.mark(boundaries, ind_zmax)

if (dimen > 2):
   ind_wave = ind_zmin
   ind_therm = ind_zmax
else:
     ind_wave = ind_bottom
     ind_therm = ind_top

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define variational problem
T  = TrialFunction(V)
vT = TestFunction(V)
T0 = Function(V)
T0 = interpolate(Constant(T_ini),V)
Ph0 = Function(V)
Ph0 = interpolate(Constant(1.0),V) 
Ph_curr = Function(V)
Ph_curr = interpolate(Constant(1.0),V)  
T_curr = Function(V)
Elmr = FiniteElement("CG", mesh.ufl_cell(), degree_L)
Elmi = FiniteElement("CG", mesh.ufl_cell(), degree_L)
Elmc = Elmr*Elmi
V2 = FunctionSpace(mesh,Elmc)
E = TrialFunction(V2)
E0 = Function(V2)
(vr, vi) = TestFunctions(V2)
Er, Ei = split(E)
Er0, Ei0 = split(E0)

f = Constant(0.0)
fileInterp = File("meshtest.pvd")
fileInterp << T0

bc_T = [DirichletBC(V, T_back, boundaries, ind_therm)]

epsR_c=[]
epsR_c.append(epsSiR_c)
epsR_c.append(epsSiR_c)
epsR_c.append(epsSiR_c)

epsI_c=[]
epsI_c.append(epsSiI_c)
epsI_c.append(epsSiI_c)
epsI_c.append(epsSiI_c)

epsR_l=[]
epsR_l.append(epsSiR_l)
epsR_l.append(epsSiR_l)
epsR_l.append(epsSiR_l)

epsI_l=[]
epsI_l.append(epsSiI_l)
epsI_l.append(epsSiI_l)
epsI_l.append(epsSiI_l)

rhoF_c=[]
rhoF_c.append(rh0Si_c)
rhoF_c.append(rh0Si_c)
rhoF_c.append(rh0Si_c)

rhoF_l=[]
rhoF_l.append(rh0Si_l)
rhoF_l.append(rh0Si_l)
rhoF_l.append(rh0Si_l)

cpF_c=[]
cpF_c.append(cpSi_c)
cpF_c.append(cpSi_c)
cpF_c.append(cpSi_c)

cpF_l=[]
cpF_l.append(cpSi_l)
cpF_l.append(cpSi_l)
cpF_l.append(cpSi_l)

kF_c=[]
kF_c.append(kSi_c)
kF_c.append(kSi_c)
kF_c.append(kSi_c)

kF_l=[]
kF_l.append(kSi_l)
kF_l.append(kSi_l)
kF_l.append(kSi_l)

LheatF_c=[]
LheatF_c.append(LHeatSi_c)
LheatF_c.append(LHeatSi_c)
LheatF_c.append(LHeatSi_c)

phequ_l=[]
phequ_l.append(1) # Si melting
phequ_l.append(0) # Si substrate
phequ_l.append(0)

# Real permittivity constant value in crystal and liquiq
def epsR_Ph(T,Ph,lay):  # Attention lay = 0 is the layer marked by 1 (0 is Air) !
    value=abs(Ph)*eval(epsR_c[lay]) + (1.0-abs(Ph))*eval(epsR_l[lay])
    return value

# Imaginary permittivity constant values in crystal and liquiq
def epsI_Ph(T,Ph,lay):
    value=abs(Ph)*eval(epsI_c[lay])+(1.0-abs(Ph))*eval(epsI_l[lay])
    return value

def cp_TAir(T):  # [J*kg-1*K-1] specific heat
    value=1047.63657-0.372589265*T+9.45304214E-4*pow(T,2)-6.02409443E-7*pow(T,3)+1.2858961E-10*pow(T,4)
    return value

def cp_Ph(T,Ph,lay):  # [J*kg-1*K-1] specific heat
    value=abs(Ph)*eval(rhoF_c[lay])*eval(cpF_c[lay]) + (1.0-abs(Ph))*eval(rhoF_l[lay])*eval(cpF_l[lay])
    return value

def k_TAir(T):
    value=-0.00227583562+1.15480022E-4*T-7.90252856E-8*pow(T,2)+4.11702505E-11*pow(T,3)-7.43864331E-15*pow(T,4)
    value=value*scale # Att.
    return value

def k_Ph(T,Ph,lay):
    slp = 0.1
    value = abs(Ph)*eval(kF_c[lay]) + (1.0-abs(Ph))*eval(kF_l[lay])
    value = value*scale
    return value

def latent(lay):  # [J / m^3] latent heat
    value=eval(rhoF_c[lay])*eval(LheatF_c[lay])
    return value

assembling = False
E_outr=Constant(1.0)
E_outi=Constant(0.0)
TE_comp = 1.0

def solver_E(Ti, 
         linear_solver='Direct',
         abs_tol=1E-5,
         rel_tol=1E-3,
         div_lim=1E4,
         no0_ini=False):
    linear_solver='Direct'
    # Define variational problem for the Maxwell equation TE
    #Q: Are those 4 elements constants during the simulation ?  #A No
    ar = inner(grad(Er), grad(vr))*dx(0) - k0sq*Er*vr*dx(0) #+ cp_Ph(expr,Ph,cpc_c,cpc_l)*Ei*vr*dx(0)
    Lr = f*vr*dx(0)
    ai = inner(grad(Ei), grad(vi))*dx(0) - k0sq*Ei*vi*dx(0) #- cp_Ph(expr,Ph,cpc_c,cpc_l)*Er*vi*dx(0)
    Li = f*vi*dx(0)
    i=1
    while i<=num_lay:
        ar = ar + inner(grad(Er), grad(vr))*dx(i) - k0sq*(epsR_Ph(Ti,Ph0,i-1)*Er*vr + epsI_Ph(Ti,Ph0,i-1)*Ei*vr)*dx(i)
        Lr = Lr + f*vr*dx(i)
        ai = ai + inner(grad(Ei), grad(vi))*dx(i) - k0sq*(epsR_Ph(Ti,Ph0,i-1)*Ei*vi - epsI_Ph(Ti,Ph0,i-1)*Er*vi)*dx(i)
        Li = Li + f*vi*dx(i)
        i+=1
    ar = ar - k0_wnumb*Ei*vr*ds(ind_wave)
    ai = ai + k0_wnumb*Er*vi*ds(ind_wave)
    Lr = Lr - 2.0*k0_wnumb*E_outi*vr*ds(ind_wave)
    Li = Li + 2.0*k0_wnumb*E_outr*vi*ds(ind_wave)
    if assembling :
        a = assemble(ar+ai)
        L = assemble(Lr+Li)
        if linear_solver == 'Krylov':
            if size_mpi == 1: # serial
                solver = KrylovSolver('gmres', 'ilu')
            else :            # parallel
                solver = KrylovSolver('gmres', 'hypre_euclid')
            solver.parameters["absolute_tolerance"] = abs_tol
            solver.parameters["relative_tolerance"] = rel_tol
            solver.parameters["divergence_limit"] = div_lim
            solver.parameters["nonzero_initial_guess"]=no0_ini
        else:
            solver = LUSolver()
        solver.solve (a, E0.vector(), L)
    else :
        problem = LinearVariationalProblem(ar+ai, Lr+Li, E0)
        solver = LinearVariationalSolver(problem)
        # Set linear solver parameters
        prm = solver.parameters
        if linear_solver == 'Krylov':
            prm["linear_solver"] = 'gmres'
            if size_mpi == 1: # serial
                prm["preconditioner"] = 'ilu'
            else :            # parellel
                prm["preconditioner"] = 'hypre_euclid'
            prm["krylov_solver"]["absolute_tolerance"] = abs_tol
            prm["krylov_solver"]["relative_tolerance"] = rel_tol
            prm["krylov_solver"]["divergence_limit"] = div_lim
            prm["krylov_solver"]["nonzero_initial_guess"]=no0_ini
        else:
            prm["linear_solver"] = 'lu'
        # Compute solution
        solver.solve()

def solverNL1D_onlyT(Ti,Phi,Phi2,dt_val,power,linear_solver='Direct',
             abs_tol=1E-5,
             rel_tol=1E-3,
             div_lim=1E4,
             no0_ini=False):
    FactS = Constant(1.) # Constant(1.e9)
    dHdt = Phi2-Phi # 30.*Phi**2*(Phi-1.0)**2*(Phi2-Phi) # Constant(0.0) # Hentalpy lost-gained left hand
    dt=Constant(dt_val)
    dtsour = Constant(dt_val*Elas*power)
    Tsol=Function(V)
    aT = cp_TAir(Ti)*T*vT*dx(0) + k_TAir(Ti)*dt*dot(grad(T), grad(vT))*dx(0)
    LT = cp_TAir(Ti)*Ti*vT*dx(0)
    i=1
    while i<=num_lay:   # Attention layer dx[0] is air while layer with index = 0 in the element list is the first layer in the structure file
        if phequ_l[i-1]==1 :
            aT = aT + cp_Ph(Ti,Phi,i-1)*T*vT*dx(i)  + dt*k_Ph(Ti,Phi,i-1)*dot(grad(T), grad(vT))*dx(i)
            LT = LT + cp_Ph(Ti,Phi,i-1)*Ti*vT*dx(i) + latent(i-1)*dHdt*vT*dx(i) + dtsour*epsI_Ph(Ti,Phi,i-1)*Trasm*vT*dx(i)
            dHdtpartial = assemble(latent(i-1)*dHdt*dx(i))
        else:
            aT = aT + cp_Ph(Ti,Phi,i-1)*T*vT*dx(i) + dt*k_Ph(Ti,Phi,i-1)*dot(grad(T), grad(vT))*dx(i)
            LT = LT + cp_Ph(Ti,Phi,i-1)*Ti*vT*dx(i) + dtsour*epsI_Ph(Ti,Phi,i-1)*Trasm*vT*dx(i)
        i+=1
    if(T_resistor) :
        aT = aT + FactS*(k_Ph(Ti,Phi,num_lay-1)/AddedSubstrate)*T*vT*ds(ind_therm)
        LT = LT + FactS*(k_Ph(Ti,Phi,num_lay-1)/AddedSubstrate)*T_back*vT*ds(ind_therm)
    if assembling :
        AT = assemble(aT)
        bT = assemble(LT)
        if(T_resistor==False) :
            for bc in bc_T : bc.apply(AT,bT)
        if linear_solver == 'Krylov':
            if size_mpi == 1: # serial
                solver = KrylovSolver('gmres', 'ilu')
            else :            # parallel
                solver = KrylovSolver('gmres', 'hypre_euclid')
            solver.parameters["absolute_tolerance"] = abs_tol
            solver.parameters["relative_tolerance"] = rel_tol
            solver.parameters["divergence_limit"] = div_lim
            solver.parameters["nonzero_initial_guess"]=no0_ini
        else:
            solver = LUSolver()
        solver.solve (AT, Tsol.vector(), bT)
    else :
        if(T_resistor==False) :
            problem = LinearVariationalProblem(aT, LT, Tsol, bc_T)
        else :
            problem = LinearVariationalProblem(aT, LT, Tsol)
        solver = LinearVariationalSolver(problem)
        # Set linear solver parameters
        prm = solver.parameters
        if linear_solver == 'Krylov':
            prm["linear_solver"] = 'gmres'
            if size_mpi == 1: # serial
                prm["preconditioner"] = 'ilu'
            else :            # parellel
                prm["preconditioner"] = 'hypre_euclid'
            prm["krylov_solver"]["absolute_tolerance"] = abs_tol
            prm["krylov_solver"]["relative_tolerance"] = rel_tol
            prm["krylov_solver"]["divergence_limit"] = div_lim
            prm["krylov_solver"]["nonzero_initial_guess"]=no0_ini
        else:
            prm["linear_solver"] = 'lu'
        # Compute solution
        solver.solve()
    return Tsol, dHdtpartial

# ----------------------------------------------------------------------------------

mpclass = Si_melting_MulskipsProcess('Si') # if not provided as argument, the calibration is internally fixed

# Box for MulSKIPS
# lenx=round((xBox1-xBox0)/alat) #nm
# leny=round((yBox1-yBox0)/alat) #nm
# lenz=round((zBox1-zBox0)/alat) #nm
lenx, leny, lenz = 480, 480, 760 # KMC superlattice units
alat = mpclass.KMC_sf *0.1 # 0.543/12. # nm ; 0.543 nm is the Silicon lattice constant
print('Box size in MulSKIPS units is:', lenx, leny, lenz)  # 480,480,1440 = (20.720, 20.720, 65.160)/alat
print('A_lat in MulSKIPS is: {} nm'.format(alat))

# Origin of Mulskips cell in fem; nm with respect to xBox0, yBox0 and zBox0
xOff, yOff, zOff = 0, 0, 196.0 # nm

# Total top air in KMC + nucleus thickness
LenVac = 40 # air thickness in KMC box (Ang)
LenNuc = 70 # air + nucleus radius (Ang)

cadfilename = 'LattGeo.dat'
tempfilename = 'TempStat.dat'
phasesfilename = 'LattStat.dat'
restartfilename = 'LattChk.dat'

# Indices for regions in geo file. 
# 0 --> air (not interpreted a wall in mulskips); 10 --> walls; 1 --> Si
#### NWtest.geo
# first region is air in geo, but we want mulskips to treat it as walls so we put 10
#regions=[10,10,1,10] # it is regions=[0,10,1,10] in geo 
#### flat.geo
regions=[10,1,10] # it is regions=[0,1,10] in geo 

# Setup geometry to be used as input in MulSKIPs; write LattGeo.dat and get cell_map, rank_map
cell_map, rank_map, wall_map = write_geom_for_mulskips(cadfilename, regions, mesh, V, dofmap, mv, comm,
                       alat, xBox0+xOff, yBox0+yOff, zBox0+zOff, lenx, leny, lenz)
# # in case you saved maps with savemaps=True and now want to save time...
# print('Reading maps...'); t2 = time.time()
# cell_map = np.loadtxt('cell_map.dat').astype(int)
# rank_map = np.loadtxt('rank_map.dat').astype(int)
# wall_map = np.loadtxt('wall_map.dat').astype(int)
# print('DONEreading maps...ETA = ', time.time() - t2)

# Melting temperature Si [K]
Tm = mpclass.Tm 

# Start multi-scale FENICS-MulSKIPs thermal simulation
execpath = '/home/gcalogero/LA/newcalib/SiKMC_nonhomo/test-merge-v2/mulskips-source-merge-v2'
count_aboveTm = 0
max_aboveTm = 1000000
# dt_val = dt_min
i_tab=0
l_tab=len(p_tab)-1
dt_val=dt_maxT
KMCstartflag = False
KMCregrowth_incomplete = True

while t <= t_stop:
    if i_tab < l_tab-1:  # interpolate power in laser file
        while ((p_tab[i_tab+1][0] < t + toffset) and (i_tab < l_tab-1)): # error while (p_tab[i_tab+1][0] < t + toffset and  i_tab < l_tab-1 ) :
            i_tab+=1
        #Q: is this still ok for QA pulse ? This kind of test probably means there is no interpolation in the pulse for QA
        #F: Linear interpolation
        if abs(p_tab[i_tab][0]-(t+toffset)) <=1.e-04:  # 1.e-13*scale
            power = p_tab[i_tab][1]
        else:
            power = (p_tab[i_tab][1]+(p_tab[i_tab+1][1]-p_tab[i_tab][1])*((t+toffset)-p_tab[i_tab][0]) \
            /(p_tab[i_tab+1][0]-p_tab[i_tab][0]))
        Source = Constant(Elas*power)
    solver_E(T0,linearsolver,abstol,reltol,divlim,no0ini)
    Er0, Ei0 = split(E0)
    Trasm=TE_comp*k0_wnumb*(Ei0*Ei0+Er0*Er0)
    T_curr, dHdtpartial=solverNL1D_onlyT(T0,Ph0,Ph_curr,dt_val,power,linearsolver,abstol,reltol,divlim,no0ini)
    T0.assign(T_curr)
    T0_nodal=T0.vector().get_local()
    T0_max=numpy.amax(T0_nodal)
    MPI.barrier(comm)
    T0_max_glob=MPI.max(comm,T0_max)
    if verb_log:
        E0_nodal=E0.vector().get_local()
        E0_max=numpy.amax(E0_nodal)
        MPI.barrier(comm)
        E0_max_glob=MPI.max(comm,E0_max)
        if rank == 0 :
            print("source evaluated TE", t, power, E0_max_glob,T0_max_glob, dHdtpartial)
        MPI.barrier(comm)

    if T0_max_glob >= Tm:      # this will activate the first time that T > Tm; nothing will restore KMCstartflag to False inside the loop, so OK
        KMCstartflag = True

    if KMCstartflag and KMCregrowth_incomplete:  
            
        dt_val=dt_min # nanoseconds
        
        # Directory for kmc calculation
        rundirname = 'KMC_t_{:.1f}'.format(t)
        
        # Save solution in VTK format
        file = File("emfield_{:.1f}.pvd".format(t))
        file << E0
        file = File("Temperature_{:.1f}.pvd".format(t))
        file << T_curr

        # Write T map for mulskips
        write_T_for_mulskips(tempfilename, t, T_curr, \
            cell_map, rank_map, mesh, V, dofmap, coord_mesh, dofs, comm, \
            alat, xBox0+xOff, yBox0+yOff, zBox0+zOff, lenx, leny, lenz)
        MPI.barrier(comm)

        # Use dt_max before turning on mulskips, then refine it 
        MaxTime_KMC = dt_val*1e-9 # dt_val in seconds
        Nframes = 100

        # OLD CALIBRATION
        # MaxTime_KMC = dt_val*1e6 # seconds
        # MaxTime_KMC = 1e5 s is equivalent to dt_min = 0.1 ns (at T>Tm) as it allows 
        # to melt roughly 0.1 nm of material. 
        # This roughly agrees with v(T) = -1 m/s given by fuchs-vogel slightly above Tm.
        # NB: Remember to fix this when correctly calibrating P0 and the energetic 
        # parameters used for probabilities (one per coor, not 0.96 only for all coor like it is now)

        # mulskips should read the original cadfilename only at the first step. 
        # Then it reads its own LattCoo
        if count_aboveTm == 0:
            # only in the first call we will set liquid in the first few top layers, 
            # otherwise mulskips does not start... 
            tmpcadname = cadfilename
            runtype = 'R'
        else:
            tmpcadname = restartfilename
            runtype = 'C' # continuation run, reading restartfilename

        # Setup and Run MulSKIPs 
        if rank==0: 
            endok = run_mulskips(execpath, rundirname, Simulation='LA', PtransZig=1.0, RunType=runtype, IDUM=9117116, \
                OutMolMol=10000, IterMax=1e15, OutTime=MaxTime_KMC/Nframes, TotTime=MaxTime_KMC, ExitStrategy='Time', mp=mpclass, \
                cadfilename=tmpcadname, tempfilename=tempfilename, phasesfilename=phasesfilename, LenVac=LenVac, LenNuc=LenNuc, \
                homogeneous=False, SaveFinalState=True, restartfilename=restartfilename)
        MPI.barrier(comm)
        
        # Check if mulskips ended with few particles error
        if not endok:
            print('\nSTOP!!! Few MC particles in KMC...probably everything has melted.')
            break

        # Read phases generated by mulskips
        Ph0 = Ph_curr
        Ph_curr, Nliquid = read_phases_from_mulskips(phasesfilename, rundirname, \
            t, cell_map, rank_map, wall_map, V, dofmap, comm, \
            alat, xBox0+xOff, yBox0+yOff, zBox0+zOff, lenx, leny, lenz) 
        MPI.barrier(comm)

        # Check if re-growth is complete
        if Nliquid < 3:
            print('\nWARNING!!! Too few liquid atoms...probably re-growth is complete.')
            KMCregrowth_incomplete = False
            dt_val=dt_maxT	
            Ph0 = Ph_curr

        count_aboveTm += 1

    t = t + dt_val

print('TOTAL ETA: {} sec'.format(time.time()-t000))

