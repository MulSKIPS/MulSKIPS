import time
import math
import numpy as np
import os,shutil,subprocess,sys

##### --------------------------------------------------------------------------
# ANALYSIS
##### --------------------------------------------------------------------------

def read_output_files(rundir, what='undercoordinated'):
    import glob

    # If what='undercoordinated': Collect all *.xyz files within the run_dir directory 
    # related to the undercoordinated atoms, that are the files I00000000.xyz, I00000001.xyz etc ...
    if what == 'undercoordinated':
        files = glob.glob(rundir+"/*[0-9].xyz")
    elif what == 'defects':
        files = glob.glob(rundir+"/*[0-9]_d.xyz")
    elif what == 'wrong':
        files = glob.glob(rundir+"/*[0-9]_w.xyz")
    elif what == 'vacancies':
        files = glob.glob(rundir+"/*[0-9]_v.xyz")
    else:
        print('ERROR:Please select one of these options for the \'what\' flag:\n \
            \'undercoordinated\', \'defects\', \'wrong\' or \'vacancies\'')
        sys.exit()

    print('There are {} \'{}\' files available in {}'.format(len(files), what, rundir))

    return files


def get_box_sides(filename):
    with open (filename) as f:
        content = f.read().splitlines()
        line=content[1]
        splitline = line.split()
        acell_x=float(splitline[1])
        acell_y=float(splitline[2])
        acell_z=float(splitline[3])
    return [acell_x,acell_y,acell_z]

def visualize(rundirname, what='undercoordinated', iteration=0, slider=False):
    """
    Visualize xyz resulting from MulSKIPs
    Note: one needs to have run the following commands to run this routine:
      jupyter nbextension install --py  widgetsnbextension --user
      jupyter nbextension enable widgetsnbextension --user --py
    """
    import ipywidgets
    import py3Dmol
    from IPython.display import display

    files = read_output_files(rundirname, what)

    if not slider:
        if iteration < len(files):
            filename = os.path.basename(files[iteration])
            print('You chose to visualize {}'.format(filename))
        else:
            print('ERROR: Iteration is too large, please set it within the range [0, {}]'.format(len(files)-1))
            sys.exit()
    else:
        print('You chose \'slider\'=True. Ignoring \'iteration\' flag.')

    def MolTo3DView(file, size=(750, 800), style="sphere", surface=False, opacity=0.5):
        """Draw molecule in 3D
        
        Args:
        ----
            xyz_file: rdMol, molecule to show
            size: tuple(int, int), canvas size
            style: str, type of drawing xyz_file
                   style can be 'line', 'stick', 'sphere', 'carton'
            surface, bool, display SAS
            opacity, float, opacity of surface, range 0.0-1.0
        Return:
        ----
            viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
        """
        assert style in ('line', 'stick', 'sphere', 'carton')
        viewer = py3Dmol.view(width=size[0], height=size[1])
        xyz_file = open(file).read()
        viewer.addModel(xyz_file,'xyz')
        if style=='sphere':
            viewer.setStyle({'sphere':{'colorscheme':'Jmol','scale':.5},'stick':{'colorscheme':'Jmol'}})
        else:
            viewer.setStyle({style:{}})
        if surface:
            viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
        viewer.rotate(-90, {'x':1,'y':0,'z':0})
        acell=get_box_sides(file)
        viewer.addLine({'start':{'x':0,'y':0,'z':0},'end':{'x':0,'y':0,'z':acell[2]}})
        viewer.addLine({'start':{'x':0,'y':0,'z':0},'end':{'x':0,'y':acell[1],'z':0}})
        viewer.addLine({'start':{'x':0,'y':0,'z':0},'end':{'x':acell[0],'y':0,'z':0}})
        viewer.addLine({'start':{'x':0,'y':acell[1],'z':acell[2]},'end':{'x':acell[0],'y':acell[1],'z':acell[2]}})
        viewer.addLine({'start':{'x':acell[0],'y':0,'z':acell[2]},'end':{'x':acell[0],'y':acell[1],'z':acell[2]}})
        viewer.addLine({'start':{'x':acell[0],'y':acell[1],'z':0},'end':{'x':acell[0],'y':acell[1],'z':acell[2]}})
        viewer.addLine({'start':{'x':0,'y':acell[1],'z':0},'end':{'x':acell[0],'y':acell[1],'z':0}})
        viewer.addLine({'start':{'x':acell[0],'y':0,'z':0},'end':{'x':acell[0],'y':acell[1],'z':0}})
        viewer.addLine({'start':{'x':0,'y':0,'z':acell[2]},'end':{'x':0,'y':acell[1],'z':acell[2]}})
        viewer.addLine({'start':{'x':0,'y':0,'z':acell[2]},'end':{'x':acell[0],'y':0,'z':acell[2]}})
        viewer.addLine({'start':{'x':0,'y':acell[1],'z':0},'end':{'x':0,'y':acell[1],'z':acell[2]}})
        viewer.addLine({'start':{'x':acell[0],'y':0,'z':0},'end':{'x':acell[0],'y':0,'z':acell[2]}})
        viewer.zoomTo()
        return viewer   

    if slider:
        def style_selector(idx, s):
            conf = files[idx]
            return MolTo3DView(conf, style=s).show()
        
        display(ipywidgets.interact(style_selector, 
                         idx=ipywidgets.IntSlider(min=0,max=len(files)-1, step=1),
                         s=ipywidgets.Dropdown(
                            options=['line', 'stick', 'sphere'],
                            value='sphere',
                            description='Style:')))
    else:
        MolTo3DView(files[iteration]).show()


# Main function of the notebook to get the surface height from a rough surface.
def get_surface_height(filename,bin_size=5.0,surface_roughness=20.0):
    """
    bin_size [Angstroem] : size of the bins in the z direction to store all z-values in an histogram
    surface_roughness [Angstroem] : surface roughness in the z direction where to average the surface height
    """
    with open(filename) as f:
        content = f.read().splitlines()

    line=content[0]
    splitline = line.split()
    nat = int(splitline[0]) # Total number of atoms
    #print('nat: {}'.format(nat))
    line=content[1]
    splitline = line.split()
    bc=str(splitline[0])
    xl=float(splitline[1])
    yl=float(splitline[2])
    zl=float(splitline[3])
    #print(zl)
    lines=content[2:]

    all_z_values = []
    all_at_species = []
    for line in lines:
        splitline = line.split()
        at=str(splitline[0])
        x=float(splitline[1])
        y=float(splitline[2])
        z=float(splitline[3])
        all_z_values.append(z)
        all_at_species.append(at)
    if len(all_z_values) != nat:
        print('Error: nat and len(all_z_values) do not coincide')
        sys.exit()

    #print(max(all_z_values),min(all_z_values))

    # if all z values are equal or too close (< binsize), do not waste time on making hyst
    if all_z_values.count(all_z_values[0]) == len(all_z_values) or (max(all_z_values)-min(all_z_values)) < bin_size: # ang
        print('* WARNING: z range in xyz file is null or < bin size, so we directly average z column.')
        return sum(all_z_values)/len(all_z_values)
    else:
        #print(zl)
        #print(zl/bin_size)
        databins=int(zl/bin_size)
        #print(databins)
        hist = np.histogram(all_z_values,bins=databins,range=(0.0, zl))
        #print(hist[0])
        maxIndexList = [i for i,j in enumerate(hist[0]) if j==max(hist[0])] #here,i=index and j = value of that index
        #print(maxIndexList)
        #print(maxIndexList[0])
        #print(hist[0][maxIndexList])
        z_min = float(bin_size * maxIndexList[0])
        z_max = float(z_min + bin_size)
        #print(z_min,z_max)
        #print(np.unique(all_z_values))
        if z_max < min(all_z_values):
            print('* WARNING: z_max < min(all_z_values), so we go to the next bin.')
            z_min+=bin_size
            z_max+=bin_size

        # average of z-value within the bin with maximum number of elements
        count = 0.0
        z_sum = 0.0
        for z in all_z_values:
            if z > z_min and z < z_max:
                count+=1
                z_sum = z_sum + z
        z_ave = float(z_sum/count)
        #print(z_ave)

        # average of z-value for atoms falling in the range (z_ave-surface_roughness,z_ave+surface_roughness)
        count = 0.0
        z_sum = 0.0
        z_min = z_ave - surface_roughness
        z_max = z_ave + surface_roughness
        for ia, z in enumerate(all_z_values):
            if z > z_min and z < z_max:
                count+=1
                z_sum = z_sum + z                     
        z_ave = float(z_sum/count)
        #print(z_ave)   

        return z_ave

def analyze_growth_rate(rundirname, bin_size=5.0,surface_roughness=20.0, method='finitediff',
    plotting=True, savepng=False, Nexclude=2, minframes=None):
    """
    Growth rate extraction
    The following notebook allows to extract the growth rate from a Super lattice 
    Kinetic Monte Carlo simulation with the mulskips code. It is supposed that you
    run mulskips for a flat (001) surface, that is "F" letter as input geometry in 
    the file start.dat (see the tutorial to run the epitaxial growth of a surface with mulskips).
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import UnivariateSpline

    # Firstly, let set the folder where you ran mulskips. 
    files = read_output_files(rundirname, 'undercoordinated')

    if minframes is None:
        minframes = Nexclude*2 +2

    # First we check if there are enough files to calculate the growth rate
    if len(files) <= minframes:
        print('ERROR: there are less than Nexclude*2 +1 files. You need more to calculate the rate.') 
        print('You should try to: i) reduce the output frequency in mulskips ii) increase TotTime or max number of iterations in mulskips iii) reduce Nexclude flag')
        print('We are now returning zero')
        return 0
    else:
        """
        To calculate the growth rate we have to extract the surface height zsurfave 
        at different KMC steps and the running time 𝑡. Then the growth rate 𝜂 is
        𝜂 = d(zsurfave)/d𝑡

        Let start with the determination of the surface height zsurfave from a 
        generic output file of the undercoordinated atoms. 
        To do that we divide the orthorhombic simulation box [x𝑙, y𝑙, z𝑙] in slices 
        along the z-direction. The size of each z-slice in angstroem is defined by 
        the variable bin_size. Essentially we are storing the z coordinates of all 
        atoms in a histogram with bins of size bin_size.

        Then we check the slice along z which contains the largest number of atoms. 
        This z-slice corresponds in a first approximation to the region where the 
        growing surface lies. We calculate the average z coordinate 𝑧largestbinave 
        for atoms belonging to such z-slice. Then we calculate zsurfave as the average 
        of the z-coordinate of all atoms falling in the range 
           |𝑧 − 𝑧largestbinave| < surface_roughness 
        To set the surface_roughness parameter it is recommended to have a look to 
        the *.xyz files by means of a visualization software to evaluate the 
        extension/quality of the surface roughness.
        """

        # All the above operation to get z surfave for a given KMC xyz file are 
        # coded in the get_surface_height function.
        # Here we calculate and store the surface height for all output files.
        all_surface_heights = []
        for file in files:
            print(file)
            surface_height  = get_surface_height(file, bin_size, surface_roughness)
            all_surface_heights.append(surface_height)
        # print(all_surface_heights)

        """
        Here we extract the KMC time from the run output. 
        When you run the mulskips code please store the screen output in the file log.txt:
        <path-of-the-compiled-code>/mulskips.e | tee log.txt
        """
        def get_time(filename):
            with open (filename) as f:
                content = f.read().splitlines()
                line=content[0]
                splitline = line.split()
                time = float(splitline[3]) # KMC time
            return time

        kmc_time_list = []
        for file in files:
            time = get_time(file)
            kmc_time_list.append(time)
        #print(kmc_time_list)

        """
        Here we have to convert the KMC time in the real process time. 
        That depends on the calibration strategy. A reasonable rescale 
        can be obtained with a constant jump frequency 𝜈=1e−12 s.
        """
        nu = 1.0 # 1.0e-12 # jump frequency
        time_list = [i * nu for i in kmc_time_list]

        # Get growth rate
        # Fit with straight line and derive
        def get_growth_rate(time_list,all_surface_heights, method='finitediff'):
            if method == 'finitediff':
                gr = []
                dx = all_surface_heights[1] - all_surface_heights[0]
                dt = time_list[1] - time_list[0]
                curr_gr = (dx / dt)*3600*1e-4 # [micron/hour]
                gr.append(curr_gr)
                for ind in range(len(time_list)-2):
                    ind = ind+1
                    dx = all_surface_heights[ind+1] - all_surface_heights[ind-1]
                    dt = time_list[ind+1] - time_list[ind-1]
                    curr_gr = (dx / dt)*3600*1e-4 # [micron/hour]
                    gr.append(curr_gr) 
                dx = all_surface_heights[-1] - all_surface_heights[-2]
                dt = (time_list[-1] - time_list[-2]) # [s]
                curr_gr = (dx / dt)*3600*1e-4 # [micron/hour] 
                gr.append(curr_gr)
                # Average
                gr_average = np.mean(np.array(gr[Nexclude:-Nexclude]))
                # gr_average = np.mean(np.array(gr[Nexclude:]))
                return gr, gr_average, None

            elif method == 'spline':
                x = np.asarray(time_list)
                y = np.asarray(all_surface_heights)
                spl = UnivariateSpline(x, y)
                der = spl.derivative()
                gr = der(x)*3600*1e-4
                fitsurfheights = spl(x)[Nexclude:-Nexclude]
                # Average
                gr_average = np.mean(np.array(gr[Nexclude:-Nexclude]))
                return gr, gr_average, fitsurfheights

            elif method == 'polyfit':
                x = np.asarray(time_list[Nexclude:-Nexclude])
                y = np.asarray(all_surface_heights[Nexclude:-Nexclude])
                c, stats = np.polynomial.polynomial.polyfit(x,y,1,full=True)
                print(stats[0][0]) # SSR should be small!
                gr_average = c[1]*3600*1e-4
                fitsurfheights = c[0] + c[1]*x
                # Average
                gr = np.ones(len(time_list))*gr_average
                return gr, gr_average, fitsurfheights

            if method == 'finiteANDfit':
                gr = []
                dx = all_surface_heights[1] - all_surface_heights[0]
                dt = time_list[1] - time_list[0]
                curr_gr = (dx / dt)*3600*1e-4 # [micron/hour]
                gr.append(curr_gr)
                for ind in range(len(time_list)-2):
                    ind = ind+1
                    dx = all_surface_heights[ind+1] - all_surface_heights[ind-1]
                    dt = time_list[ind+1] - time_list[ind-1]
                    curr_gr = (dx / dt)*3600*1e-4 # [micron/hour]
                    gr.append(curr_gr) 
                dx = all_surface_heights[-1] - all_surface_heights[-2]
                dt = (time_list[-1] - time_list[-2]) # [s]
                curr_gr = (dx / dt)*3600*1e-4 # [micron/hour] 
                gr.append(curr_gr)
                # Average
                gr_average = np.mean(np.array(gr[Nexclude:-Nexclude]))
                # gr_average = np.mean(np.array(gr[Nexclude:]))
                
                # Fit
                x = np.asarray(time_list[Nexclude:-Nexclude])
                y = np.asarray(all_surface_heights[Nexclude:-Nexclude])
                c, stats = np.polynomial.polynomial.polyfit(x,y,1,full=True)
                print(stats[0][0]) # SSR should be small!
                gr_average_fit = c[1]*3600*1e-4
                fitsurfheights = c[0] + c[1]*x

                return gr, (gr_average, gr_average_fit), fitsurfheights

        growth_rate, gr_ave, mysurfheights = get_growth_rate(time_list,all_surface_heights,method)

        print('\nAverage growth rate for Process ID {}: {} [micron/hour]\n'.format(rundirname,gr_ave))
        # print('Experimental growth rates')
        # for key in data_samples:
        #     print('Process ID {}: {} [micron/hour]'.format(key,data_samples[key]['Exp-Growth-Rate']))

        # Plot stuff
        if plotting:
            # Here we plot the surface height as a zsurfave  as function of the KMC steps and the process time
            plt.figure(figsize=(10,4))
            plt.rcParams.update({'font.size': 14})
            plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            plt.subplot(121)
            plt.plot(time_list,all_surface_heights, 'k-')
            if method != 'finitediff':
                plt.plot(time_list[Nexclude:-Nexclude],mysurfheights, 'r--')
            plt.ylabel('Surface height [$\AA$]')
            plt.xlabel('time [s]')
            # plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

            plt.subplot(122)
            plt.plot(time_list,growth_rate, 'k-')
            plt.plot(time_list[Nexclude:-Nexclude],growth_rate[Nexclude:-Nexclude], 'r-')
            # plt.plot(time_list[Nexclude:],growth_rate[Nexclude:], 'r-')
            plt.ylabel('Growth rate $\eta$ [micron/hour]')
            plt.xlabel('time [s]')
            # plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            
            # plt.title('Process ID: {}'.format(rundirname[:5]))        
            plt.tight_layout()
            if savepng:
                #rundir = os.getcwd() +'/'+rundirname+'/'
                #plt.savefig('fig-{}-surface_height-growth_rate-time.png'.format(rundirname[:5]))
                plt.savefig('fig-{}.png'.format(rundirname))
            plt.show()

        return gr_ave
            

# Main function of the notebook to get the coverage out of the xyz
def get_coverage(filename, mp=None):
    """
    CAREFUL: here only H is considered. It Needs some adjustments for whene other cov species are included 
    """
    with open(filename) as f:
        content = f.read().splitlines()

    line=content[0]
    splitline = line.split()
    nat = int(splitline[0]) # Total number of atoms
    #print('nat: {}'.format(nat))
    line=content[1]
    splitline = line.split()
    bc=str(splitline[0])
    xl=float(splitline[1])
    yl=float(splitline[2])
    zl=float(splitline[3])
    #print(zl)
    lines=content[2:]

    all_at_species = []
    for line in lines:
        splitline = line.split()
        at=str(splitline[0])
        x=float(splitline[1])
        y=float(splitline[2])
        z=float(splitline[3])
        all_at_species.append(at)
    if len(all_at_species) != nat:
        print('Error: nat and len(all_at_species) do not coincide')
        sys.exit()
    
    coverage_H, coverage_Cl, nvac, nundercoo = 0, 0, 0, 0
    for ia in all_at_species:
        # count number of coverages around surface
        if ia == 'H':
            coverage_H += 1
        elif ia == 'Cl':
            coverage_Cl += 1
        elif ia == 'O':
            nvac += 1
        else:  # this is all undercoordinated crystal atoms (including those surrounding vacancies)
            nundercoo += 1

    # We do not consider as "surface" the 4 crystal atoms surrounding each vacancy (O atoms in xyz) 
    surface_atoms = nundercoo - nvac*4   

    # Get number of dangling bonds per surface atom
    if mp.miller == 100:
        at2dang = 1.0 #but remember that the max theoretical coverage for Si(100) (num of DBs per Si) is 2; it could be changed in case it's needed
    else:
        at2dang = 1.0 #max theoretical coverage for Si(110) and (111) (num of DBs per Si)
    
    # Get number of dangling bonds in the KMC box
    surface_dang_bonds = surface_atoms * at2dang
    
    # Get coverage (it should be between 0 and 1)
    coverage = {'H': [coverage_H / surface_dang_bonds], 'Cl': [coverage_Cl / surface_dang_bonds]}

    # Instead this is how it would be if normalized to the number of surface dangling bonds in a FLAT surface
    # area_box_KMC = xl*yl*1e-20   # m2
    # coverage = coverage_H / (mp.at_density * area_box_KMC * at2dang)
    
    return coverage 



def analyze_coverage(rundirname, plotting=True, savepng=False, Nexclude=2, minframes=None, mp=None):
    
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import UnivariateSpline

    if mp is None:
        print('ERROR: please provide process class argument (mp) in analyze_growth_rate...')
        sys.exit()

    # Firstly, let set the folder where you ran mulskips. 
    files = read_output_files(rundirname, 'undercoordinated')

    # First we check if there are enough files to calculate the growth rate
    if minframes is None:
        minframes = Nexclude*2 +2
    if len(files) <= minframes:
        print('ERROR: there are less than Nexclude*2 +1 files. You need more to calculate the rate.') 
        print('You should try to: i) reduce the output frequency in mulskips ii) increase TotTime or max number of iterations in mulskips iii) reduce Nexclude flag')
        print('We are now returning zero')
        return 0
    else:
        # All the above operation to get z surfave for a given KMC xyz file are 
        # coded in the get_surface_height function.
        # Here we calculate and store the surface height for all output files.
        # all_surface_coverage = {}
        for i, file in enumerate(files):
            print(file)
            surface_coverage = get_coverage(file, mp)
            if i==0: 
                all_surface_coverage = surface_coverage.copy() # will be zeros because first file is I00000000.xyz
            else:
                for key in surface_coverage:
                    all_surface_coverage[key].extend(surface_coverage[key])
        #print(all_surface_coverage)

        """
        Here we extract the KMC time from the run output. 
        When you run the mulskips code please store the screen output in the file log.txt:
        <path-of-the-compiled-code>/mulskips.e | tee log.txt
        """
        def get_time(filename):
            with open (filename) as f:
                content = f.read().splitlines()
                line=content[0]
                splitline = line.split()
                time = float(splitline[3]) # KMC time
            return time

        kmc_time_list = []
        for file in files:
            time = get_time(file)
            kmc_time_list.append(time)
        #print(kmc_time_list)

        """
        Here we have to convert the KMC time in the real process time. 
        That depends on the calibration strategy. A reasonable rescale 
        can be obtained with a constant jump frequency 𝜈=1e−12 s.
        """
        nu = 1.0 # 1.0e-12 # jump frequency
        time_list = [i * nu for i in kmc_time_list]

        """
        1) conto il numero N_h di idrogeni in un'intorno (+- roughness) della average surface height in ogni xyz di un run mulskips
        2) il numero di dangling bonds dei silici in superficie, assumendo flat surface, è N_Si100 = area_box_KMC * 6.78e14 at/cm2
        3) coverage = N_h / N_Si100
        """
        coverage_ave = {}
        for key in all_surface_coverage:
            coverage_ave[key] = np.mean(np.asarray(all_surface_coverage[key][Nexclude:-Nexclude]))
            print('\nAverage {} coverage for Process ID {}: {} \n'.format(key, rundirname, coverage_ave[key]))

        # Plot stuff
        if plotting:
            colors = {'Cl':'b', 'H':'g'}
            # Here we plot the surface height as a zsurfave  as function of the KMC steps and the process time
            plt.figure(figsize=(10,4))
            plt.rcParams.update({'font.size': 14})
            plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            for key in all_surface_coverage:
                plt.plot(time_list, all_surface_coverage[key], '-o', color=colors[key], label=key, alpha=0.5)
                plt.axhline(coverage_ave[key], color=colors[key], ls=':')
            plt.ylabel('Coverage')
            plt.xlabel('time [s]')
            
            # plt.title('Process ID: {}'.format(rundirname[:5]))        
            plt.legend(loc=0)
            plt.tight_layout()
            if savepng:
                #rundir = os.getcwd() +'/'+rundirname+'/'
                #plt.savefig('fig-{}-surface_height-growth_rate-time.png'.format(rundirname[:5]))
                plt.savefig('fig-{}.png'.format(rundirname))
            plt.show()

        return coverage_ave


def export_xyz(xyzfile, newfile, alat, what='surface+coverage', DEP3Dfile=None, mesh=None, voff=[0,0,0], reverse_z=False, exclude=[]):
    
    if not what in ['surface', 'coverage', 'surface+coverage']:
        print("ERROR: what should one of the following: \n'surface', 'coverage', 'surface+coverage'")
        sys.exit() 

    if DEP3Dfile is not None:
        if mesh is None:
            print('ERROR: missing argument for Dolfin mesh in export_xyz()')
            sys.exit()

    # Read MulSKIPS output
    spec, x, y, z = np.loadtxt(xyzfile, skiprows=2, unpack=True,
        dtype={'names':('spec', 'x', 'y', 'z'), 'formats':('U2', 'f8', 'f8', 'f8')})
    xyz = np.array(list(zip(x,y,z)))
    print('Number of atoms in MulSKIPS output: ', len(xyz))

    """
    Small filtering to fix internal bug in MulSKIPS.
    Sometimes coverage atoms are duplicates of vacancies (O in file).
    The lines below ensure that the duplicate coverages are removed.
    The indices returned by np.unique will retain only O 
    NB: we reverse the list to have all O at the beginning of the list.
    This ensures that the second found duplicate is always the wrong coverage.
    """
    spec, xyz = spec[::-1], xyz[::-1] # this is enough to put ox at the beginning of the list
    a, idx = np.unique(xyz, axis=0, return_index=True)
    spec, xyz = spec[idx], xyz[idx]
    print('Number of atoms in MulSKIPS output (after removing duplicates): ', len(xyz))

    # Indices separated per species 
    Hlist = (spec == 'H').nonzero()[0]
    Cllist = (spec == 'Cl').nonzero()[0]
    Olist = (spec == 'O').nonzero()[0]
    Silist = (spec == 'Si').nonzero()[0]
    Gelist = (spec == 'Ge').nonzero()[0]

    if what == 'coverage': 
    # Write only coordinates of coverage atoms (coincides with surface for high coverages)
        xyz_final = xyz[np.concatenate((Hlist, Cllist))]
        spec_final = spec[np.concatenate((Hlist, Cllist))]
    
    else:  # Coordinates of coverage + surface or only surface (see below) 
        from scipy import spatial
        bondlength = alat * (3**0.5)/4 # Ang
        toremove = []

        # First filter the system without translation
        tree = spatial.KDTree(xyz)
        for i in Olist:
            dd, iout = tree.query(xyz[i], k=5) # O atom itself + its max 4 neighbours
            toremove += iout[dd < bondlength+0.1].tolist()

        # Still Some Si atoms wrapped due to PBC are not captured 
        # So now repeat the strategy after translating a bit the structure along x and y
        # Remember to wrap all atoms falling outside back inside the box
        xyz1 = xyz + [2*bondlength, 2*bondlength, 0]
        xmax, ymax = np.amax(xyz[:,0]), np.amax(xyz[:,1])
        xyz1[xyz1[:,0] > xmax] -= [xmax, 0, 0] # wrap back inside along x
        xyz1[xyz1[:,1] > ymax] -= [0, ymax, 0] # wrap back inside along y
        tree = spatial.KDTree(xyz1)
        for i in Olist:
            dd, iout = tree.query(xyz1[i], k=5) # O atom itself + its max 4 neighbours
            toremove += iout[dd < bondlength+0.1].tolist()

        if what == 'surface':
            toremove += Hlist.tolist()
            toremove += Cllist.tolist()

        # Now remove all of them (no need to remove duplicates in toremove list) 
        xyz_final = np.delete(xyz, toremove, axis=0)
        spec_final = np.delete(spec, toremove, axis=0)

    print(f'Number of {what} atoms: ', len(xyz_final))

    # Exclude atoms if needed
    if len(exclude) > 0: # zero-based
        xyz_final = np.delete(xyz_final, exclude, axis=0)
        spec_final = np.delete(spec_final, exclude, axis=0)

    # Write output XYZ file just to visualize
    with open(xyzfile) as fin:
        head = [next(fin) for x in range(2)]
    np.savetxt(newfile, np.c_[spec_final, xyz_final], fmt='%s', 
        header=str(len(spec_final))+"\n"+head[-1].strip('\n'), comments='')

    # Write output XYZ file ready to be used in DEP3D (micron units)
    if DEP3Dfile is not None:

        if reverse_z: # if LA then I need to reverse z, because MulKSIPS has surface facing positive z, while mesh is the other way around 
            xyz_final = xyz_final*[1,1,-1]

        xyz = mesh.coordinates() # nm
        for ii in range(3):
            xyz_final[:,ii] += (xyz[:,ii].min() + voff[ii]) * 10 # Ang
            
        print('Writing DEP3D file after shifting coordinates onto initial mesh')
        np.savetxt(DEP3Dfile, xyz_final*1e-4, fmt='%.8f', comments='', header=str(len(xyz_final)))

        # import the xyz below in paraview together with *_KMCregions* pvd file (both in nm!) 
        cell = np.array(head[-1].strip('\n').split()[1:]).astype(np.float32) # Ang
        np.savetxt(newfile+"aligned.xyz", np.c_[spec_final, (xyz_final*0.1).astype(np.float32)], fmt='%s', 
            header=str(len(spec_final))+"\nperiodic "+' '.join((cell*0.1).astype(str)), comments='')

    return


