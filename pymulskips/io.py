from dolfin import *
import sys, time
import numpy as np
import os
from subprocess import call
import meshio
import cashocs  

def msh2dolfin(mshname, save_h5=True, rescale=None, rotate_angle=None, rotate_axis=None):
	
	structurename = mshname[:-4]

	# Print this to check what's inside msh
	print('\nINFO msh file:\n', meshio.read(f"{structurename}.msh")) 

	# Remove existing xdmf files
	try:
	    os.remove(f"{structurename}_boundaries.xdmf")
	except OSError:
	    pass
	try:
	    os.remove(f"{structurename}_subdomains.xdmf")
	except OSError:
	    pass
	# Convert msh to xdmf and h5. This will produce 3+3 files (or 2+2 if there are no boundaries)
	call([f"cashocs-convert {structurename}.msh {structurename}.xdmf"], shell=True)

	# Import xdmf and h5 to FEniCS
	mesh, subdomains, boundaries, dx, ds, dS = cashocs.import_mesh(f'{structurename}.xdmf')

	# Rescale/rotate mesh if needed
	### WE ASSUME TESTMUNDFAB.MSH IS IN mm !!!! Check it with Remi...or get a smaller one...
	if rescale is not None:
	    mesh.scale(rescale)
	if rotate_angle != None and rotate_axis != None:
	    mesh.rotate(rotate_angle, rotate_axis)
	    
	# Read mesh
	if save_h5:
		hdf = HDF5File(mesh.mpi_comm(), f"{structurename}_all.h5", "w")
		hdf.write(mesh, "/mesh")
		hdf.write(subdomains, "/subdomains")
		hdf.write(boundaries, "/boundaries")
		hdf.close()

	return mesh, subdomains, boundaries, dx, ds, dS


def dolfin2mulskips(cadfilename, mpclass, regions, mesh, subdomains, 
    regionsfilename='KMCregions', return_maps=True, savemaps=False):

    comm = MPI.comm_world
    if(MPI.rank(comm)==0):
        print('\nStarted dolfin2mulskips...')
        t0 = time.time()
    MPI.barrier(comm)

    # Check regions
    print('Setting the following regions in MulSKIPS:')
    for reg in regions.items():
        print(f'Subdomain {reg[0]} in mesh \t --> Region of type {reg[1]} in MulSKIPS')

    # KMC lattice parameter in nm
    alat = mpclass.KMC_sf * 1e-1

    # Box for MulSKIPS should be determined directly from the mesh size
    # Mulskips should then be compiled accordingly
    xyz = mesh.coordinates() # nm
    Lx = np.abs(xyz[:,0].max() - xyz[:,0].min())
    Ly = np.abs(xyz[:,1].max() - xyz[:,1].min())
    Lz = np.abs(xyz[:,2].max() - xyz[:,2].min())
    print('alat in MulSKIPS = {} nm'.format(alat))
    print('Box size in msh (nm) =', Lx, Ly, Lz)  
    # KMC superlattice units alat (nm)
    lenx = int(round(Lx/alat))
    leny = int(round(Ly/alat))
    lenz = int(round(Lz/alat)) 
    print('Box size in MulSKIPS units (using Si alat=0.04525 nm) =', lenx, leny, lenz)
    lenx -= lenx%24
    leny -= leny%24
    print('Box size to be set in MulSKIPS to ensure periodicity =', lenx, leny, lenz)
    print('      NB: MulSKIPS box strain vector w.r.t. the original mesh: ', (lenx*alat - Lx)/Lx, (leny*alat - Ly)/Ly, (lenz*alat - Lz)/Lz )

    xBox0 = xyz[:,0].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges 
    yBox0 = xyz[:,1].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges
    zBox0 = xyz[:,2].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges
    print('xBox0, yBox0, zBox0 = ', xBox0, yBox0, zBox0)

    # Get mesh and dofmap 
    degree_L=1
    lagr_el = FiniteElement('CG', mesh.ufl_cell(), degree_L)
    V = FunctionSpace(mesh, lagr_el) #c
    dofmap = V.dofmap()

    lattcell_map = np.zeros((lenz,leny*lenx), dtype=int)
    cell_map = np.zeros((lenz,leny*lenx), dtype=int)
    rank_map = np.ones((lenz,leny*lenx), dtype=int)*(MPI.size(comm)+1)  # to avoid confusion with rank numbers

    # Write to file to check that regions are correctly set up 
    my_first, my_last = dofmap.ownership_range() # global
    LattFunc = Function(V)
    LattFunc.assign(Constant(0))
    LattFunc_array = LattFunc.vector().get_local()
    for cell in cells(mesh):
        ind_cur = subdomains[cell] # -1
        if(regions[ind_cur] == 10):
          dofs_cell = dofmap.cell_dofs(cell.index())  # local
          for dof in dofs_cell :
              global_dof = dofmap.local_to_global_index(dof)  # global
              if my_first <= global_dof < my_last:
                 LattFunc_array[dof]=10
    for cell in cells(mesh):
        ind_cur = subdomains[cell]
        if(regions[ind_cur] == 1):
          dofs_cell = dofmap.cell_dofs(cell.index())  # local
          for dof in dofs_cell :
              global_dof = dofmap.local_to_global_index(dof)  # global
              if my_first <= global_dof < my_last:
                 LattFunc_array[dof]=1
    MPI.barrier(comm)
    LattFunc.vector().set_local(LattFunc_array)
    LattFunc.vector().apply('insert')
    file = File(regionsfilename+".pvd")
    file << LattFunc

    # Tranform information on regions from mesh to mulskips superlattice
    if(MPI.rank(comm)==0): print('Tranform information on regions from mesh to mulskips superlattice...')
    MPI.barrier(comm)
    bbox_tree, ncells = mesh.bounding_box_tree(), mesh.num_cells()
    
    # cell_map_flat = cell_map.ravel()
    # rank_map_flat = rank_map.ravel()
    # lattcell_map_flat = lattcell_map.ravel()
    # t1 = time.time()
    # for izj in range(lenz*leny*lenx):
    #     if(izj%(10*leny*lenx) == 0 and MPI.rank(comm)==0):
    #         print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(int(izj/(lenx*leny)), lenz, time.time()-t1), end='')
    #         t1 = time.time()
    #     # 1d index: izj = ix + iy*lenx + iz*lenx*leny
    #     ix = izj % (lenx-1)
    #     iy = (izj / (lenx-1)) % (leny-1)
    #     iz = izj / ((lenx-1)*(leny-1))
    #     point_latt=Point(xyz[:,0].min()+alat*float(ix), xyz[:,1].min()+alat*float(iy), xyz[:,2].min()+alat*float(iz))
    #     val = -1
    #     cell_index = bbox_tree.compute_first_entity_collision(point_latt) # it gives infty if out of thread subdomain
    #     if cell_index < ncells:
    #         val = regions[subdomains[cell_index]]
    #         cell_map_flat[izj] = cell_index
    #         rank_map_flat[izj] = MPI.rank(comm)
    #     MPI.barrier(comm)
    #     val = int(MPI.max(comm, val))
    #     lattcell_map_flat[izj] = val
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)
    # cell_map = cell_map_flat.reshape((lenz,leny*lenx))
    # rank_map = rank_map_flat.reshape((lenz,leny*lenx))
    # lattcell_map = lattcell_map_flat.reshape((lenz,leny*lenx))

    # the lines above are slightly faster but give a slightly different LattGeo.dat 
    iz = 0
    # LattFunc.set_allow_extrapolation(True)
    if(MPI.rank(comm)==0): t1 = time.time()
    while iz < lenz :
      iy = 0
      if(iz%10 == 0 and MPI.rank(comm)==0):
        print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(iz, lenz, time.time()-t1), end='')
        t1 = time.time()
      while iy < leny :
        ix = 0
        while ix < lenx :
          j=ix+iy*lenx
          point_latt=Point(xBox0+alat*float(ix), yBox0+alat*float(iy), zBox0+alat*float(iz))
          val = -1
          cell_index = bbox_tree.compute_first_entity_collision(point_latt) # it gives infty if out of thread subdomain
          if cell_index < ncells:
            val = regions[subdomains[cell_index]]
            cell_map[iz][j] = cell_index
            rank_map[iz][j] = MPI.rank(comm)
          MPI.barrier(comm)
          val = int(MPI.max(comm, val))
          if val == -1:
            print(ix, iy, iz, point_latt.x(), point_latt.y(), point_latt.z(), cell_index, val)
            print('Stopping here...'); sys.exit() # uncomment to check
    #      val = LattFunc(point_latt)
          lattcell_map[iz][j] = val
          ix+=1
        # if iy==2:
        #     print('Stopping here...'); sys.exit() # uncomment to check
        iy+=1
      iz+=1
    if(MPI.rank(comm)==0):print('')
    MPI.barrier(comm)

    if(MPI.rank(comm)==0): print('DONE Setting up geometry for mulskips. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    # Write LattGeo.dat input for MulSKIPs
    if(MPI.rank(comm)==0):
        print('Writing geometry to {}...'.format(cadfilename))
        t0 = time.time()
        fileout=open(cadfilename,"w")
        lines=0
        col=lenz-1
        while col >= 0:
            while lines < leny*lenx:
                fileout.write("%s " % str(lattcell_map[col][lines]))
                lines+=1
            fileout.write("\n")
            col-=1
            lines=0
        fileout.close()
        print('DONE Writing geometry to {}. ETA: {} sec'.format(cadfilename, time.time()-t0))
    MPI.barrier(comm)

    if savemaps:
        print('Saving cell_map, rank_map and wall_map to file for restart...'); t0 = time.time()
        np.savetxt('cell_map_{}.dat'.format(MPI.rank(comm)), cell_map, fmt='%d')
        np.savetxt('rank_map_{}.dat'.format(MPI.rank(comm)), rank_map, fmt='%d')
        np.savetxt('wall_map_{}.dat'.format(MPI.rank(comm)), LattFunc_array, fmt='%d')
        print('DONE Saving maps files. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    if return_maps:
        return cell_map, rank_map, LattFunc_array


def mulskips2dolfin(coofilename, mpclass, mesh, subdomains, regions,
    cell_map, rank_map, wall_map, LA=False, 
    pvdfilename='Mesh_AfterKMC'):
    
    comm = MPI.comm_world
    if(MPI.rank(comm)==0):
        print('\nStarted mulskips2dolfin...')
        t0 = time.time()
    MPI.barrier(comm)

    # # Check regions
    # print('Setting the following regions in MulSKIPS:')
    # for reg in regions.items():
    #     print(f'Subdomain {reg[0]} in mesh \t --> Region of type {reg[1]} in MulSKIPS')

    # KMC lattice parameter in nm
    alat = mpclass.KMC_sf * 1e-1

    # Box for MulSKIPS should be determined directly from the mesh size
    # Mulskips should then be compiled accordingly
    xyz = mesh.coordinates() # nm
    Lx = np.abs(xyz[:,0].max() - xyz[:,0].min())
    Ly = np.abs(xyz[:,1].max() - xyz[:,1].min())
    Lz = np.abs(xyz[:,2].max() - xyz[:,2].min())
    print('alat in MulSKIPS = {} nm'.format(alat))
    print('Box size in msh (nm) =', Lx, Ly, Lz)  
    # KMC superlattice units alat (nm)
    lenx = int(round(Lx/alat))
    leny = int(round(Ly/alat))
    lenz = int(round(Lz/alat)) 
    print('Box size in MulSKIPS units =', lenx, leny, lenz)
    lenx -= lenx%24
    leny -= leny%24
    print('Box size to be set in MulSKIPS to ensure periodicity =', lenx, leny, lenz)

    xBox0 = xyz[:,0].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges 
    yBox0 = xyz[:,1].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges
    zBox0 = xyz[:,2].min() + 1e-10 # this 1e-10 is to avoid numerical imprecisions when colliding the point in the mesh at edges
    print('xBox0, yBox0, zBox0 = ', xBox0, yBox0, zBox0)

    # Read coofilename input from MulSKIPs
    print('Reading Coor (LA phases) from {}...'.format(coofilename))
    t0 = time.time()
    value_latt  = np.zeros((lenz,leny*lenx), dtype=int)

    with open(coofilename,"r") as filein: # faster than np.loadtxt
        lines = filein.readlines()
        for i, line in enumerate(reversed(lines)):
            words = line.split()
            value_latt[i] = words
    MPI.barrier(comm)
    if len(value_latt[0]) != lenx*leny: # Check
        print('ERROR: Dimension of coofilename is not ({}, {})'.format(lenz,leny*lenx))
        sys.exit()
    print('DONE Reading Coor (LA phases) from {}. ETA: {} sec'.format(coofilename, time.time()-t0))

    # #-----------------------------------
    # # DEBUG: this should be identical to LattStat.dat!!! 
    # print('Writing geometry TESTTTTTTTTT to {}...'.format('test.dat'))
    # t0 = time.time()
    # fileout=open('test.dat',"w")
    # lines=0
    # col=lenz-1
    # while col >= 0:
    #     while lines < leny*lenx:
    #         fileout.write("%s " % str(value_latt[col][lines]))
    #         lines+=1
    #     fileout.write("\n")
    #     col-=1
    #     lines=0
    # fileout.close()
    # print('DONE Writing geometry TESTTTTTTTTT to {}. ETA: {} sec'.format('test.dat', time.time()-t0)) 
    # #------------------------------------

    if(MPI.rank(comm)==0): print('Mapping sites from MulSKIPs superlattice to fem mesh...')
    MPI.barrier(comm)
    # Loop over cells, get relative dofs and cumulate value reading from superlattice
    # Note: dofs can belong to multiple cells.
    t0 = time.time()

    # Get mesh and dofmap 
    degree_L=1
    lagr_el = FiniteElement('CG', mesh.ufl_cell(), degree_L)
    V = FunctionSpace(mesh, lagr_el) #c
    dofmap = V.dofmap()

    # Initialize new fem vector of zeros 
    Ph_tmp = Function(V)
    Ph_tmp.assign(Constant(0.0))
    Ph_tmp_array = Ph_tmp.vector().get_local()

    # Find unique cell indices in cell_map which correspond to the KMC box 
    cell_map_flat = cell_map.ravel()
    KMCcells, iKMCc = np.unique(cell_map_flat, return_inverse=True)
    # Note that KMCcells[iKMCc] == cell_map_flat --> print(np.allclose(KMCcells[iKMCc], cell_map_flat))
    # We will use this below to map back the dofs

    # Now we find the map of dofs corresponding to the unique(cell_map_flat)!
    KMC_dofs_cell = [] # it will always be sized [len(KMCcells), 4] if not MPI. 
    # List is probably better than array when MPI is used, and when dof number x cell varies..
    # For speed gain, initialize directly as np.zeros(len(KMCcells), 4)
    my_first, my_last = dofmap.ownership_range() # global
    for ic, kmccell in enumerate(KMCcells):   # this loop is muuuch faster than those over KMC sites
        # print('\r', 'Iter {:7d} /{:7d}'.format(ic, len(KMCcells)), end='')
        dofs_cell = dofmap.cell_dofs(kmccell)  # local
        dofs_cell_ok = []
        for dof in dofs_cell:
            global_dof = dofmap.local_to_global_index(dof)  # global
            if my_first <= global_dof < my_last:
                dofs_cell_ok.append(dof)
        KMC_dofs_cell.append(dofs_cell_ok)
    MPI.barrier(comm)
    # if(MPI.rank(comm)==0):print(''); MPI.barrier(comm)   <-------- QUESTO è IL MALE. NON FARLO MAI.

    # Now we just need to expand the map back to the lenz*lenx*leny size, using the iKMCc indices
    dofs_cell_map = np.asarray(KMC_dofs_cell)[iKMCc]
    # This array now contains a 1:1 correspondence between each KMC site and its corresponding dofs in the mesh
    # It can be used to select elements from Ph_tmp_array!

    # We now need to cumulate value_latt elements on the corresponding dofs within Ph_tmp_array 
    # We will restrict only to where value_latt is different from 10 or 0 
    # (we already know where walls are from wall_map, and we don't care about zeros because we are only modifying Ph_tmp_array where it is solid)
    value_latt_flat = value_latt.ravel()
    KMCnotwall = np.isin(value_latt_flat, [1,2,3,4]).nonzero()[0] 

    # For 1 processor, this is how we can fill Ph_tmp_array in one line: 
    # Ph_tmp_array[dofs_cell_map[KMCnotwall]] += value_latt_flat[KMCnotwall].reshape(-1, 1)  # reshape is needed just to allow broadcasting
    # As an alternative, the block below works also for MPI, and it is fast!
    # ----------------------------------------
    rank_map_flat = rank_map.ravel()
    for izj in KMCnotwall:
        if MPI.rank(comm) == rank_map_flat[izj]:
            Ph_tmp_array[dofs_cell_map[izj]] += value_latt_flat[izj]
        MPI.barrier(comm)
    MPI.barrier(comm)
    # ----------------------------------------

    # Set =0 if the dof is liquid or wall, =1 if dof is solid
    Ph_tmp_array[Ph_tmp_array > 0] = 1.  # solid=1, rest is 0 for now...(remember that Ph_tmp_array does not contain 10!)
    nsolid = len((Ph_tmp_array > 0).nonzero()[0])

    if not LA:
        Ph_tmp_array[wall_map > 5] = 10.   # now walls are = 10, and only liquid/air is =0
    else:
        Ph_tmp_array[wall_map > 5] = 1.   # now walls are also =1, and only liquid is =0
    nliquid = len((Ph_tmp_array == 0).nonzero()[0])

    # Print out info
    nPh = len(Ph_tmp_array)
    print('Process {}: Total: {}'.format(MPI.rank(comm), nPh))
    MPI.barrier(comm)
    print('Process {}: Crystal sites: {}'.format(MPI.rank(comm), nsolid))
    MPI.barrier(comm)
    print('Process {}: Air / Liquid sites: {}'.format(MPI.rank(comm), nliquid))
    MPI.barrier(comm)
    print('Process {}: Non-evolving sites: {}'.format(MPI.rank(comm), nPh-nsolid-nliquid))

    # Sum over processors
    MPI.barrier(comm)
    n_tot = int(MPI.sum(comm, nPh))
    nsolid_tot = int(MPI.sum(comm, nsolid))
    nliquid_tot = int(MPI.sum(comm, nliquid))
    if MPI.rank(comm) == 0:
        print('Total Crystal sites: {}'.format(nsolid_tot))
        print('Total Air / Liquid sites: {}'.format(nliquid_tot))
        print('Total Non-evolving sites: {}'.format(n_tot-nsolid_tot-nliquid_tot))
    MPI.barrier(comm)

    # # Savetxt if you need to double check...
    # np.savetxt('Ph_tmp_array_p{}.dat'.format(MPI.rank(comm)), Ph_tmp_array)
    # print('Process {} has saved Ph_tmp_array of type {} to file... '.format(MPI.rank(comm), type(Ph_tmp_array)))
    # MPI.barrier(comm)
    

    # # TO BENCHMARK COMPARE TO THIS ----------------------------------------------------------------------------
    # t0 = time.time()
    # # Initialize new fem vector of zeros 
    # Ph_tmp_2 = Function(V)
    # Ph_tmp_2 = interpolate(Constant(0.0), V)
    # Ph_tmp_array_2 = Ph_tmp_2.vector().get_local()
    # # Loop over cells, get relative dofs and cumulate value reading from superlattice
    # # Note: dofs can belong to multiple cells.
        
    # # faster than 3 nested loops
    # cell_map_flat = cell_map.ravel()
    # rank_map_flat = rank_map.ravel()
    # value_latt_flat = value_latt.ravel()
    # t1 = time.time()
    # for izj in range(lenz*leny*lenx):
    #     if(izj%(10*leny*lenx) == 0 and MPI.rank(comm)==0):
    #         print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(int(izj/(lenx*leny)), lenz, time.time()-t1), end='')
    #         t1 = time.time()
    #     if value_latt_flat[izj] != 10: 
    #         if MPI.rank(comm) == rank_map_flat[izj]:
    #             dofs_cell = dofmap.cell_dofs(cell_map_flat[izj])  # local
    #             Ph_tmp_array_2[dofs_cell] += value_latt_flat[izj]
    #         MPI.barrier(comm)
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)

    # # Set =0 if the dof is liquid or wall, =1 if dof is solid
    # Ph_tmp_array_2[Ph_tmp_array_2 > 0] = 1.  # solid=1, rest is 0 for now...
    # nsolid = len((Ph_tmp_array_2 > 0).nonzero()[0])
    # Ph_tmp_array_2[wall_map > 5] = 1.   # now walls are also =1, and only liquid is =0

    # Ph_tmp_2.vector().set_local(Ph_tmp_array_2)
    # Ph_tmp_2.vector().apply('insert')
    # print('DONE Mapping phases from MulSKIPs superlattice to fem mesh. ETA: {} sec'.format(time.time()-t0))

    # # Check how many are solid 
    # nliquid = len((Ph_tmp_array_2 == 0).nonzero()[0])
    # if MPI.rank(comm) == 0:
    #     print('Solid phase atoms: {}'.format(nsolid))
    #     print('Liquid phase atoms: {}'.format(nliquid))
    #     print('Walls: {}'.format(len(Ph_tmp_array_2)-nsolid-nliquid))
    # MPI.barrier(comm)

    # print('Is Ph_tmp_array the same as before optimization?', np.allclose(Ph_tmp_array, Ph_tmp_array_2))
    # print((Ph_tmp_array == Ph_tmp_array_2).sum(), '/', len(Ph_tmp_array), 'are the same!')
    # # TO BENCHMARK COMPARE TO THIS ----------------------------------------------------------------------------

    if(MPI.rank(comm)==0): print('DONE Mapping MulSKIPs superlattice to fem mesh. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    if not LA:
        # Map back from KMC regions to MSH subdomain indices
        mapped_values = np.ones(len(Ph_tmp_array), dtype=int)*1000
        for mshID, kmcID in regions.items():
            mapped_values[Ph_tmp_array == kmcID] = mshID

        # Overwrite subdomains indices
        subdomains_new_values = np.ones(mesh.num_cells(), dtype=int)*1000
        for cell in cells(mesh):   # this loop is muuuch faster than those over KMC sites
            icell = cell.index()
            dofs_cell = dofmap.cell_dofs(icell)  # local
            for dof in dofs_cell:
                global_dof = dofmap.local_to_global_index(dof)  # global
                if my_first <= global_dof < my_last:
                    subdomains_new_values[icell] = mapped_values[dof]
        subdomains.set_values(subdomains_new_values)

    # FirOverwrite original meshfunction 
    Ph_tmp.vector().set_local(Ph_tmp_array)
    Ph_tmp.vector().apply('insert')

    ### Visualize only solid and walls as 1, otherwise 0
    if(MPI.rank(comm)==0): print(f'Writing phases {pvdfilename}.pvd file...'); t0 = time.time()
    MPI.barrier(comm)
    fileInterp = File(pvdfilename+".pvd")
    fileInterp << Ph_tmp
    if(MPI.rank(comm)==0): print(f'DONE writing phases {pvdfilename}.pvd file. ETA: {time.time()-t0} sec')
    MPI.barrier(comm)

    ### Visualize final subdomains
    if(MPI.rank(comm)==0): print(f'Writing phases {pvdfilename}_subdomains.pvd file...'); t0 = time.time()
    MPI.barrier(comm)
    fileInterp = File(pvdfilename+"_subdomains.pvd")
    fileInterp << subdomains
    if(MPI.rank(comm)==0): print(f'DONE writing phases {pvdfilename}_subdomains.pvd file. ETA: {time.time()-t0} sec')
    MPI.barrier(comm)

    if not LA:
        return subdomains
    else:
        return Ph_tmp, nliquid_tot


### Seems that it works, but can be refined, there must be a better way of writing it...
def dolfin2msh(in_mshname, out_mshname, mesh, subdomains, boundaries):

    # Remember to convert subdomains elements into normal int32 (they are uint64!)
    subarray = subdomains.array().astype(int)
    print(subarray)
    
    # Read input msh
    M = meshio.read(in_mshname)

    # NO NEED TO ROTATE OR RESCALE BECAUSE THE MESH WILL BE READ FROM THE INITIAL ONE!
    def create_mesh(mesh, cell_type, new_cell_data):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        # Double Check
        if new_cell_data.shape != cell_data.shape:
            print('new_cell_data.shape != cell_data.shape.Aborting...')
            sys.exit()
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, \
            cell_data={"Data_From_MulSKIPS":[new_cell_data]})
        return out_mesh
    tetra_mesh = create_mesh(M, "tetra", new_cell_data=subarray)
    tetra_mesh.write('tmp.msh', binary=False)

    # Appending ElementData of out_msh to in_msh
    from collections import deque
    with open(in_mshname) as fp:
        data = fp.read()
    with open('tmp.msh') as fp:
        data2 = ''.join(list(deque(fp,len(subarray)+10))) # this is like UNIX tail
    data += data2  
    with open (out_mshname, 'w') as fp:
        fp.write(data)














# ---------------------- LA routines -------------------------

def write_geom_for_mulskips(cadfilename, regions, mesh, V, dofmap, mv, comm,
    alat, xBox0, yBox0, zBox0, lenx, leny, lenz, savemaps=False):

    if(MPI.rank(comm)==0):
        print('Setting up geometry for mulskips...')
        t0 = time.time()
    MPI.barrier(comm)

    lattcell_map = np.zeros((lenz,leny*lenx), dtype=int)
    cell_map = np.zeros((lenz,leny*lenx), dtype=int)
    rank_map = np.ones((lenz,leny*lenx), dtype=int)*(MPI.size(comm)+1)  # to avoid confusion with rank numbers

    # Write to file to check that regions are correctly set up 
    my_first, my_last = dofmap.ownership_range() # global
    LattFunc = Function(V)
    LattFunc.assign(Constant(0))
    LattFunc_array = LattFunc.vector().get_local()
    for cell in cells(mesh):
        ind_cur = mv[cell]
        if(regions[ind_cur] == 10):
          dofs_cell = dofmap.cell_dofs(cell.index())  # local
          for dof in dofs_cell :
              global_dof = dofmap.local_to_global_index(dof)  # global
              if my_first <= global_dof < my_last:
                 LattFunc_array[dof]=10
    for cell in cells(mesh):
        ind_cur = mv[cell]
        if(regions[ind_cur] == 1):
          dofs_cell = dofmap.cell_dofs(cell.index())  # local
          for dof in dofs_cell :
              global_dof = dofmap.local_to_global_index(dof)  # global
              if my_first <= global_dof < my_last:
                 LattFunc_array[dof]=1
    MPI.barrier(comm)
    LattFunc.vector().set_local(LattFunc_array)
    LattFunc.vector().apply('insert')
    file = File("LattFunc.pvd")
    file << LattFunc

    # Tranform information on regions from mesh to mulskips superlattice
    if(MPI.rank(comm)==0): print('Tranform information on regions from mesh to mulskips superlattice...')
    MPI.barrier(comm)
    bbox_tree, ncells = mesh.bounding_box_tree(), mesh.num_cells()
    
    # cell_map_flat = cell_map.ravel()
    # rank_map_flat = rank_map.ravel()
    # lattcell_map_flat = lattcell_map.ravel()
    # t1 = time.time()
    # for izj in range(lenz*leny*lenx):
    #     if(izj%(10*leny*lenx) == 0 and MPI.rank(comm)==0):
    #         print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(int(izj/(lenx*leny)), lenz, time.time()-t1), end='')
    #         t1 = time.time()
    #     # 1d index: izj = ix + iy*lenx + iz*lenx*leny
    #     ix = izj % (lenx-1)
    #     iy = (izj / (lenx-1)) % (leny-1)
    #     iz = izj / ((lenx-1)*(leny-1))
    #     point_latt=Point(xBox0+alat*float(ix), yBox0+alat*float(iy), zBox0+alat*float(iz))
    #     val = -1
    #     cell_index = bbox_tree.compute_first_entity_collision(point_latt) # it gives infty if out of thread subdomain
    #     if cell_index < ncells:
    #         val = regions[mv[cell_index]]
    #         cell_map_flat[izj] = cell_index
    #         rank_map_flat[izj] = MPI.rank(comm)
    #     MPI.barrier(comm)
    #     val = int(MPI.max(comm, val))
    #     lattcell_map_flat[izj] = val
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)
    # cell_map = cell_map_flat.reshape((lenz,leny*lenx))
    # rank_map = rank_map_flat.reshape((lenz,leny*lenx))
    # lattcell_map = lattcell_map_flat.reshape((lenz,leny*lenx))

    # the lines above are slightly faster but give a slightly different LattGeo.dat 
    iz = 0
    # LattFunc.set_allow_extrapolation(True)
    if(MPI.rank(comm)==0): t1 = time.time()
    while iz < lenz :
      iy = 0
      if(iz%10 == 0 and MPI.rank(comm)==0):
        print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(iz, lenz, time.time()-t1), end='')
        t1 = time.time()
      while iy < leny :
        ix = 0
        while ix < lenx :
          j=ix+iy*lenx
          point_latt=Point(xBox0+alat*float(ix), yBox0+alat*float(iy), zBox0+alat*float(iz))
          val = -1
          cell_index = bbox_tree.compute_first_entity_collision(point_latt) # it gives infty if out of thread subdomain
          if cell_index < ncells:
            val = regions[mv[cell_index]]
            cell_map[iz][j] = cell_index
            rank_map[iz][j] = MPI.rank(comm)
          MPI.barrier(comm)
          val = int(MPI.max(comm, val))
    #      val = LattFunc(point_latt)
          lattcell_map[iz][j] = val
          ix+=1
        iy+=1
      iz+=1
    if(MPI.rank(comm)==0):print('')
    MPI.barrier(comm)

    if(MPI.rank(comm)==0): print('DONE Setting up geometry for mulskips. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    # Write LattGeo.dat input for MulSKIPs
    if(MPI.rank(comm)==0):
        print('Writing geometry to {}...'.format(cadfilename))
        t0 = time.time()
        fileout=open(cadfilename,"w")
        lines=0
        col=lenz-1
        while col >= 0:
            while lines < leny*lenx:
                fileout.write("%s " % str(lattcell_map[col][lines]))
                lines+=1
            fileout.write("\n")
            col-=1
            lines=0
        fileout.close()
        print('DONE Writing geometry to {}. ETA: {} sec'.format(cadfilename, time.time()-t0))
    MPI.barrier(comm)

    if savemaps:
        print('Saving cell_map, rank_map and wall_map to file for restart...'); t0 = time.time()
        np.savetxt('cell_map_{}.dat'.format(MPI.rank(comm)), cell_map, fmt='%d')
        np.savetxt('rank_map_{}.dat'.format(MPI.rank(comm)), rank_map, fmt='%d')
        np.savetxt('wall_map_{}.dat'.format(MPI.rank(comm)), LattFunc_array, fmt='%d')
        print('DONE Saving maps files. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    return cell_map, rank_map, LattFunc_array


def write_T_for_mulskips(tempfilename, tstep, T_curr, 
    cell_map, rank_map, mesh, V, dofmap, coord_mesh, dofs, comm,
    alat, xBox0, yBox0, zBox0, lenx, leny, lenz):

    # Map temperature on MulSKIPS superlattice
    if(MPI.rank(comm)==0): print('Mapping temperature on MulSKIPs superlattice...')
    MPI.barrier(comm)

    # Average T within every cell in the mesh 
    t0 = time.time()
    T_curr_array = T_curr.vector().get_local()
    print('rank/num_cells', MPI.rank(comm), mesh.num_cells())
    # if (MPI.rank(comm)==0):
    #     print('Sum num_cells', MPI.sum(comm, mesh.num_cells()))
    # MPI.barrier(comm)
    
    T_cell_av = np.zeros(mesh.num_cells())
    # T_cell_av = np.zeros(MPI.sum(comm, mesh.num_cells()))
    my_first, my_last = dofmap.ownership_range() # global
    for cell in cells(mesh):   # this loop is muuuch faster than those over KMC sites
        icell = cell.index()
        ncount = 0
        dofs_cell = dofmap.cell_dofs(icell)  # local
        for dof in dofs_cell:
            global_dof = dofmap.local_to_global_index(dof)  # global
            if my_first <= global_dof < my_last:
                T_cell_av[icell] += T_curr_array[dof]
                ncount += 1
        T_cell_av[icell] /= ncount
    MPI.barrier(comm)

    # Now the mapping from mesh to KMC lattice can be (instantly!) done with this:
    cell_map_flat = cell_map.ravel()
    rank_map_flat = rank_map.ravel()

    # This works for single core
    # value_latt = T_cell_av[cell_map_flat[rank_map_flat == MPI.rank(comm)]] # shape = (lenz*lenx*leny)

    # This might work for multi core...because value_latt keeps the shape lenz*leny*lenx and is filled by different processors at different slices  
    value_latt = np.zeros(lenz*leny*lenx)
    value_latt[rank_map_flat == MPI.rank(comm)] = T_cell_av[cell_map_flat[rank_map_flat == MPI.rank(comm)]] # shape = (lenz*lenx*leny)
    MPI.barrier(comm)

    value_latt = value_latt.reshape((lenz,leny*lenx))  # shape = (lenz, lenx*leny)
    # NB: even though T_cell_av and cell_map have different shapes and sizes, 
    # numpy will directly reshape T_cell_av correctly

    if(MPI.rank(comm)==0): print('DONE Mapping temperature on MulSKIPs superlattice. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)

    # # ----------------------------------------------------------------------------------------------
    # # Alternative un-optimized method 2:
    # t0=time.time()
    # value_latt_2  = np.zeros((lenz,leny*lenx))
    # t1 = time.time()
    # iz = 0
    # while iz < lenz:
    #   iy = 0
    #   if(iz%10 == 0 and MPI.rank(comm)==0):
    #     # print('\r', 'Iter {:7d} /{:7d}'.format(iz, lenz), end='')
    #     print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(iz, lenz, time.time()-t1), end='')
    #     t1 = time.time()
    #   while iy < leny:
    #     ix = 0
    #     while ix < lenx:
    #       j=ix+iy*lenx
    #       # if lattcell_map[iz][j] != 10: # this could speed up a lot...would it be ok?
    #       if MPI.rank(comm) == rank_map[iz][j]:
    #         dofs_cell = dofmap.cell_dofs(cell_map[iz][j])  # local
    #         value_latt_2[iz][j] = np.mean(T_curr_array[dofs_cell])
    #         # value_latt[iz][j] = T_cell_av[cell_map[iz][j]]
    #       MPI.barrier(comm)
    #       ix+=1
    #     iy+=1
    #   iz+=1
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)
    # print('DONE Mapping temperature on MulSKIPs superlattice. ETA: {} sec'.format(time.time()-t0))

    # # ----------------------------------------------------------------------------------------------
    # # Alternative un-optimized method 3: 
    # # faster than 3 nested loops, but waaaaay slower than method 1, 
    # # because T averaging is uselessly repeated many times!
    # t0 = time.time()
    # cell_map_flat = cell_map.ravel()
    # rank_map_flat = rank_map.ravel()
    # value_latt_3 = np.zeros(lenz*leny*lenx)
    # t1 = time.time()
    # for izj in range(lenz*leny*lenx):
    #     if(izj%(10*leny*lenx) == 0 and MPI.rank(comm)==0):
    #         print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(int(izj/(lenx*leny)), lenz, time.time()-t1), end='')
    #         t1 = time.time()
    #     if MPI.rank(comm) == rank_map_flat[izj]:
    #         dofs_cell = dofmap.cell_dofs(cell_map_flat[izj])  # local
    #         value_latt_3[izj] = sum(T_curr_array[dofs_cell])/len(dofs_cell)
    #         # value_latt[izj] = T_cell_av[cell_map_flat[izj]]   
    #     MPI.barrier(comm)
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)
    # value_latt_3 = value_latt_3.reshape((lenz,leny*lenx))
    # print('DONE Mapping temperature on MulSKIPs superlattice. ETA: {} sec'.format(time.time()-t0))

    # print('Are methods 1 and 2 equal?', np.array_equal(value_latt, value_latt_2))
    # print('Are methods 1 and 3 equal?', np.array_equal(value_latt, value_latt_3))
    
    # Write to file
    if(MPI.rank(comm)==0):
        print('Writing temperature map to {}...'.format(tempfilename))
        t0 = time.time()
        fileout=open(tempfilename,"w")
        lines, col = 0, lenz-1
        while col >= 0:
            while lines < leny*lenx:
                fileout.write("%s " % str(int(round(value_latt[col][lines]))))
                lines+=1
            fileout.write("\n")
            col-=1
            lines=0
        fileout.close()
        print('DONE Writing temperature to {}. ETA: {} sec'.format(tempfilename, time.time()-t0))
    MPI.barrier(comm)

    # ################  Debugging  ##########################    # NOT WORKING WHEN KMC IS ON A SUBBOX
    # print('Writing TestTemperature_{:.1f}.pvd file...'.format(tstep))
    # t0 = time.time()
    # C1 = Function(V)
    # C1 = interpolate(Constant(0.0), V)
    # C1_array = C1.vector().get_local()
    # for dof in dofs :
    #     point=coord_mesh[dof]
    #     icx=int((point[0]-xBox0)/alat)
    #     icy=int((point[1]-yBox0)/alat)
    #     icz=int((point[2]-zBox0)/alat)
    #     iz=icz
    #     if(iz<0): iz=iz+lenz
    #     if(iz>lenz-1): iz=iz-lenz
    #     iy=icy
    #     if(iy<0):iy=iy+leny
    #     if(iy>leny-1):iy=iy-leny
    #     ix=icx                       
    #     if(ix<0):ix=ix+lenx
    #     if(ix>lenx-1):ix=ix-lenx    
    #     j=ix+iy*lenx
    #     C1_array[dof] = value_latt[iz][j]
          
    # MPI.barrier(comm)    
    # C1.vector().set_local(C1_array)
    # C1.vector().apply('insert')
    # fileInterp = File("TestTemperature_{:.1f}.pvd".format(tstep))
    # fileInterp << C1
    # print('DONE Writing TestTemperature_{:.1f}.pvd file. ETA: {} sec'.format(tstep, time.time()-t0))

    return


def read_phases_from_mulskips(coofilename, rundirname, tstep, 
    cell_map, rank_map, wall_map, V, dofmap, comm,
    alat, xBox0, yBox0, zBox0, lenx, leny, lenz):
    
    # Read coofilename input from MulSKIPs
    print('Reading phases from {}...'.format(coofilename))
    t0 = time.time()
    value_latt  = np.zeros((lenz,leny*lenx), dtype=int)

    with open(coofilename,"r") as filein: # faster than np.loadtxt
        lines = filein.readlines()
        for i, line in enumerate(reversed(lines)):
            words = line.split()
            value_latt[i] = words
    MPI.barrier(comm)
    if len(value_latt[0]) != lenx*leny: # Check
        print('ERROR: Dimension of coofilename is not ({}, {})'.format(lenz,leny*lenx))
        sys.exit()
    print('DONE Reading phases from {}. ETA: {} sec'.format(coofilename, time.time()-t0))

    # #-----------------------------------
    # # DEBUG: this should be identical to LattStat.dat!!! 
    # print('Writing geometry TESTTTTTTTTT to {}...'.format('test.dat'))
    # t0 = time.time()
    # fileout=open('test.dat',"w")
    # lines=0
    # col=lenz-1
    # while col >= 0:
    #     while lines < leny*lenx:
    #         fileout.write("%s " % str(value_latt[col][lines]))
    #         lines+=1
    #     fileout.write("\n")
    #     col-=1
    #     lines=0
    # fileout.close()
    # print('DONE Writing geometry TESTTTTTTTTT to {}. ETA: {} sec'.format('test.dat', time.time()-t0)) 
    # #------------------------------------

    if(MPI.rank(comm)==0): print('Mapping phases from MulSKIPs superlattice to fem mesh...')
    MPI.barrier(comm)
    # Loop over cells, get relative dofs and cumulate value reading from superlattice
    # Note: dofs can belong to multiple cells.
    t0 = time.time()
    # Initialize new fem vector of zeros 
    Ph_tmp = Function(V)
    Ph_tmp.assign(Constant(0.0))
    Ph_tmp_array = Ph_tmp.vector().get_local()

    # Find unique cell indices in cell_map which correspond to the KMC box 
    cell_map_flat = cell_map.ravel()
    KMCcells, iKMCc = np.unique(cell_map_flat, return_inverse=True)
    # Note that KMCcells[iKMCc] == cell_map_flat --> print(np.allclose(KMCcells[iKMCc], cell_map_flat))
    # We will use this below to map back the dofs

    # Now we find the map of dofs corresponding to the unique(cell_map_flat)!
    KMC_dofs_cell = [] # it will always be sized [len(KMCcells), 4] if not MPI. 
    # List is probably better than array when MPI is used, and when dof number x cell varies..
    # For speed gain, initialize directly as np.zeros(len(KMCcells), 4)
    my_first, my_last = dofmap.ownership_range() # global
    for ic, kmccell in enumerate(KMCcells):   # this loop is muuuch faster than those over KMC sites
        # print('\r', 'Iter {:7d} /{:7d}'.format(ic, len(KMCcells)), end='')
        dofs_cell = dofmap.cell_dofs(kmccell)  # local
        dofs_cell_ok = []
        for dof in dofs_cell:
            global_dof = dofmap.local_to_global_index(dof)  # global
            if my_first <= global_dof < my_last:
                dofs_cell_ok.append(dof)
        KMC_dofs_cell.append(dofs_cell_ok)
    MPI.barrier(comm)
    # if(MPI.rank(comm)==0):print(''); MPI.barrier(comm)   <-------- QUESTO è IL MALE. NON FARLO MAI.

    # Now we just need to expand the map back to the lenz*lenx*leny size, using the iKMCc indices
    dofs_cell_map = np.asarray(KMC_dofs_cell)[iKMCc]
    # This array now contains a 1:1 correspondence between each KMC site and its corresponding dofs in the mesh
    # It can be used to select elements from Ph_tmp_array!

    # We now need to cumulate value_latt elements on the corresponding dofs within Ph_tmp_array 
    # We will restrict only to where value_latt is different from 10 or 0 
    # (we already know where walls are from wall_map, and we don't care about zeros because we are only modifying Ph_tmp_array where it is solid)
    value_latt_flat = value_latt.ravel()
    KMCnotwall = np.isin(value_latt_flat, [1,2,3,4]).nonzero()[0] 

    # For 1 processor, this is how we can fill Ph_tmp_array in one line: 
    # Ph_tmp_array[dofs_cell_map[KMCnotwall]] += value_latt_flat[KMCnotwall].reshape(-1, 1)  # reshape is needed just to allow broadcasting
    # As an alternative, the block below works also for MPI, and it is fast!
    # ----------------------------------------
    rank_map_flat = rank_map.ravel()
    for izj in KMCnotwall:
        if MPI.rank(comm) == rank_map_flat[izj]:
            Ph_tmp_array[dofs_cell_map[izj]] += value_latt_flat[izj]
        MPI.barrier(comm)
    MPI.barrier(comm)
    # ----------------------------------------

    # Set =0 if the dof is liquid or wall, =1 if dof is solid
    Ph_tmp_array[Ph_tmp_array > 0] = 1.  # solid=1, rest is 0 for now...
    nsolid = len((Ph_tmp_array > 0).nonzero()[0])
    Ph_tmp_array[wall_map > 5] = 1.   # now walls are also =1, and only liquid is =0
    nliquid = len((Ph_tmp_array == 0).nonzero()[0])

    # Print out info
    nPh = len(Ph_tmp_array)
    print('Process {}: Total: {}'.format(MPI.rank(comm), nPh))
    MPI.barrier(comm)
    print('Process {}: Solid phase atoms: {}'.format(MPI.rank(comm), nsolid))
    MPI.barrier(comm)
    print('Process {}: Liquid phase atoms: {}'.format(MPI.rank(comm), nliquid))
    MPI.barrier(comm)
    print('Process {}: Walls: {}'.format(MPI.rank(comm), nPh-nsolid-nliquid))

    # Sum over processors
    MPI.barrier(comm)
    n_tot = int(MPI.sum(comm, nPh))
    nsolid_tot = int(MPI.sum(comm, nsolid))
    nliquid_tot = int(MPI.sum(comm, nliquid))
    if MPI.rank(comm) == 0:
        print('Total Solid phase atoms: {}'.format(nsolid_tot))
        print('Total Liquid phase atoms: {}'.format(nliquid_tot))
        print('Total Walls: {}'.format(n_tot-nsolid_tot-nliquid_tot))
    MPI.barrier(comm)

    # # Savetxt if you need to double check...
    # np.savetxt('Ph_tmp_array_p{}.dat'.format(MPI.rank(comm)), Ph_tmp_array)
    # print('Process {} has saved Ph_tmp_array of type {} to file... '.format(MPI.rank(comm), type(Ph_tmp_array)))
    # MPI.barrier(comm)
    
    Ph_tmp.vector().set_local(Ph_tmp_array)
    Ph_tmp.vector().apply('insert')
    
    if(MPI.rank(comm)==0): print('DONE Mapping phases from MulSKIPs superlattice to fem mesh. ETA: {} sec'.format(time.time()-t0))
    MPI.barrier(comm)


    # # TO BENCHMARK COMPARE TO THIS ----------------------------------------------------------------------------
    # t0 = time.time()
    # # Initialize new fem vector of zeros 
    # Ph_tmp_2 = Function(V)
    # Ph_tmp_2 = interpolate(Constant(0.0), V)
    # Ph_tmp_array_2 = Ph_tmp_2.vector().get_local()
    # # Loop over cells, get relative dofs and cumulate value reading from superlattice
    # # Note: dofs can belong to multiple cells.
        
    # # faster than 3 nested loops
    # cell_map_flat = cell_map.ravel()
    # rank_map_flat = rank_map.ravel()
    # value_latt_flat = value_latt.ravel()
    # t1 = time.time()
    # for izj in range(lenz*leny*lenx):
    #     if(izj%(10*leny*lenx) == 0 and MPI.rank(comm)==0):
    #         print('\r', 'Iter {:7d} /{:7d} \t ETA x 10: {:.3f} s'.format(int(izj/(lenx*leny)), lenz, time.time()-t1), end='')
    #         t1 = time.time()
    #     if value_latt_flat[izj] != 10: 
    #         if MPI.rank(comm) == rank_map_flat[izj]:
    #             dofs_cell = dofmap.cell_dofs(cell_map_flat[izj])  # local
    #             Ph_tmp_array_2[dofs_cell] += value_latt_flat[izj]
    #         MPI.barrier(comm)
    # if(MPI.rank(comm)==0):print('')
    # MPI.barrier(comm)

    # # Set =0 if the dof is liquid or wall, =1 if dof is solid
    # Ph_tmp_array_2[Ph_tmp_array_2 > 0] = 1.  # solid=1, rest is 0 for now...
    # nsolid = len((Ph_tmp_array_2 > 0).nonzero()[0])
    # Ph_tmp_array_2[wall_map > 5] = 1.   # now walls are also =1, and only liquid is =0

    # Ph_tmp_2.vector().set_local(Ph_tmp_array_2)
    # Ph_tmp_2.vector().apply('insert')
    # print('DONE Mapping phases from MulSKIPs superlattice to fem mesh. ETA: {} sec'.format(time.time()-t0))

    # # Check how many are solid 
    # nliquid = len((Ph_tmp_array_2 == 0).nonzero()[0])
    # if MPI.rank(comm) == 0:
    #     print('Solid phase atoms: {}'.format(nsolid))
    #     print('Liquid phase atoms: {}'.format(nliquid))
    #     print('Walls: {}'.format(len(Ph_tmp_array_2)-nsolid-nliquid))
    # MPI.barrier(comm)

    # print('Is Ph_tmp_array the same as before optimization?', np.allclose(Ph_tmp_array, Ph_tmp_array_2))
    # print((Ph_tmp_array == Ph_tmp_array_2).sum(), '/', len(Ph_tmp_array), 'are the same!')
    # # TO BENCHMARK COMPARE TO THIS ----------------------------------------------------------------------------


    ### Visualize only solid as 1, otherwise 0
    if(MPI.rank(comm)==0): print('Writing phases {}/InterpPhase_{:.1f}.pvd file...'.format(rundirname, tstep)); t0 = time.time()
    MPI.barrier(comm)
    fileInterp = File(rundirname+"/InterpPhase_{:.1f}.pvd".format(tstep))
    fileInterp << Ph_tmp
    if(MPI.rank(comm)==0): print('DONE writing phases {}/InterpPhase_{:.1f}.pvd file. ETA: {} sec'.format(rundirname, tstep, time.time()-t0))
    MPI.barrier(comm)


    return Ph_tmp, nliquid_tot

# ----------------------------------------------------------------------------------

