import time
import math
import numpy as np
import os,shutil,subprocess,sys

##### --------------------------------------------------------------------------
# SETUP AND RUN MULSKIPS
##### --------------------------------------------------------------------------

# Function to Generate input for mulskips and run mulskips simulation 
def run_mulskips(execpath, runpath, Simulation, mp, PtransZig, 
    ExitStrategy='Iter', OutMolMol=None, IterMax=None, TotTime=None, OutTime=None,
    Seed_box=[120,120,120], RunType='R', IDUM=9117116, setup_only=False,
    cadfilename=None, tempfilename=None, LenVac=None, LenNuc=None, homogeneous=False,
    SaveCoo=False, coofilename=None,
    SaveFinalState=False, restartfilename=None):
    """
    Generate input for mulskips and run mulskips simulation.
    """
    t0 = time.time()

    # detect the current working directory and print it
    path = os.getenv("PWD")
    print ("The current working directory is %s" % path)


    print('\n********************************** Generating MulSKIPs input file (start.dat) **********************************')
    print('*********** Input arguments:')
    print('Executable directory: ', execpath)
    print('Simulation directory: ', runpath)
    necessary_input = [execpath, runpath, Simulation, mp, PtransZig]
    if not all(v is not None for v in necessary_input):
        print('\nERROR. Please define all the necessary input (see list below):')
        print('[execpath, runpath, Simulation, mp, PtransZig]')
        print(necessary_input)
        sys.exit()    
    print('Mulskips class: ', mp)
    print('Type of simulation: ', Simulation)
    print('PtransZig: ', PtransZig)
    print("ExitStrategy: ", ExitStrategy)
    if ExitStrategy == 'Iter':
        necessary_input = [OutMolMol, IterMax]
        if not all(v is not None for v in necessary_input):
            print('\nERROR. Please define all the necessary input (see list below):')
            print('[OutMolMol, IterMax]')
            print(necessary_input)
            sys.exit()    
        OutMolMol = int(round(OutMolMol))
        IterMax = int(round(IterMax))
        print('OutMolMol: ', OutMolMol)
        print('IterMax: ', IterMax)
    if ExitStrategy == 'Time':
        necessary_input = [OutTime, TotTime]
        if not all(v is not None for v in necessary_input):
            print('\nERROR. Please define all the necessary input (see list below):')
            print('[OutTime, TotTime]')
            print(necessary_input)
            sys.exit()    
        print('TotTime: ', TotTime, 'sec')
        print('OutTime: ', OutTime, 'sec')
    if Simulation in ['LA', 'IN']:
        print('Size of the system: from CAD')
    else:
        print('Size of the system: ', Seed_box)
    print('RunType: ', RunType)
    print('IDUM: ', IDUM)
    print('SaveFinalState: ', SaveFinalState)
    print('SaveCoo: ', SaveCoo)
    print('Only generate start.dat and do not run: ', setup_only)

    if SaveFinalState:
        savefinalstateflag = 'T'
        if restartfilename is None:
            print('ERROR: Please provide restartfilename... ')
            sys.exit()
    else:
        savefinalstateflag = 'F'

    if SaveCoo:
        savecooflag = 'T'
        if coofilename is None:
            print('ERROR: Please provide coofilename... ')
            sys.exit()
    else:
        savecooflag = 'F'

    if RunType == 'C':
        if restartfilename is None:
            print('ERROR: Please provide restartfilename... ')
            sys.exit()

    # ----- GENERATE INPUT FILE

    if mp.__class__.__name__ == 'PhaseChange':
        
        # If Laser annealing
        if Simulation == 'LA':
            """
            Generate input for mulskips and run simulation with T map and geometry from fem.
            NB: We'll set only energetics in start.dat. Probabilities will be calculated 
            directly inside mulskips using T map from fem.
            """
            # Check input is OK
            necessary_input = [TotTime, LenVac, LenNuc, tempfilename, cadfilename, SaveCoo, coofilename, restartfilename]
            if not all(v is not None for v in necessary_input):
                print('\nERROR. Please define all the necessary input (see list below):')
                print('[TotTime, LenVac, LenNuc, tempfilename, cadfilename, SaveCoo, coofilename, restartfilename]')
                print(necessary_input)
                sys.exit()
            print('Top air thickness : {} Ang'.format(LenVac))
            print('Nucleus radius : {} Ang'.format(LenNuc))
            print('Input CAD geometry file: ', path+'/'+cadfilename)    
            print('Input temperature map file: ', path+'/'+tempfilename)    
            if not SaveCoo:
                print('ERROR: Please set SevaCoo to \'T\' and provide coofilename... ')
                sys.exit()
            print('Output Coor file: ', path+'/'+coofilename)    
            print('Output checkpoint file for restart run: ', path+'/'+restartfilename)    
            if homogeneous:
                homoflag = 'T'
            else:
                homoflag = 'F'
            print('Model homogeneous nucleation: ', homoflag)

            # Generation of the evaporation energetics
            # The probabilities will be calculated internally in mulskips as follows:
            # PtransE_Si = P0 * np.exp(-ESi/(kb*Temp))
            Eevap = mp.get_energetics_evap()

            # Generation of deposition energetics
            # The probabilities will be calculated internally in mulskips as follows:
            # PtransD = 0.5 * P0 * np.exp(-Edep/(kb*Tm)) * dampingfactor(temp)
            Edep = mp.get_energetics_depo()

            # For BAck compatibility:
            PtransE = Eevap
            PtransD = Edep

        # Si(100) melting/solidification
        elif Simulation == 'FS': #or FG or FSG
            """
            Generate input for mulskips and run simulation with T map and geometry from fem.
            NB: We'll set only energetics in start.dat. Probabilities will be calculated 
            directly inside mulskips using T map from fem.
            """

            # EVAPORATION
            PtransE = mp.get_PtransE()
            # DEPOSITION
            PtransD = mp.get_PtransD()

    # in the future I will dop the samw check as below also for SL and LA
    elif mp.__class__.__name__ == 'CVD':
        # Get deposition probabilities
        PtransD = mp.get_PtransD()
        print('Deposition probabilities: \n {}'.format(PtransD))

        # Get evaporation probabilities
        """ The logic is as follows:
        We set the 3 evaporation energies for 1coor, 2coor, and 3coor Si atoms starting from the cases without coverage neighbours. 
        These are set by perturbing the energetics for deposition, with the constraint that energy for evaporation of Si 2coor is = to energy 
        for deposition of Si 1coor. The evap energies for 1coor and 3coor Si atoms are then properly calibrated.
        These 3 energies, in turn, are perturbated by the presence of coverage neighbours, by an amount, ab-initio or calibrated, 
        set by the function deltaEcov.
        Then, to find the probabilities, we use as prefactors the same fluxes used to get the deposition frequencies.
        """
        PtransE = mp.get_PtransE()
        print('Evaporation probabilities: \n {}'.format(PtransE))

        # Get absorption probabilities
        PtransAbs = mp.get_PtransAbs()
        print('Absorption probabilities: \n {}'.format(PtransAbs))

        # Get desorption probabilities
        PtransDes = mp.get_PtransDes()
        print('Desorption probabilities: \n {}'.format(PtransDes))

    elif mp.__class__.__name__ == 'PVD':
        """
        The stuff below should belong to PVD class as get_ene, get_PtransE etc, same as for CVD...
        """

        # Get evaporation probabilities
        PtransE = mp.get_PtransE()
        print('Evaporation probabilities: \n {}'.format(PtransE))

        # Get deposition probabilities
        PtransD = mp.get_PtransD()
        print('Deposition probabilities: \n {}'.format(PtransD))


    # ---------------------------------
    # Generation of the kinetic Monte Carlo input file
    def print_file_input():
        file = open('start.dat',"w") 
        NCrystal, NCov = len(mp.listcry), len(mp.listcov)
        file.write('{} {} ! NCrystal, NCov \n'.format(NCrystal, NCov))
        file.write('{} ! Atomic numbers of crystal species \n'.format(' '.join([str(x) for x in mp.listcryZ])))
        if NCov!=0:
            file.write('{} ! Atomic numbers of coverage species \n'.format(' '.join([str(x) for x in mp.listcovZ])))
        file.write(Simulation+" ! Initstat: S Sphere, C Parellelepipid, F Flat SiC-3C(100) surface, FS Flat Si(100) surface, A APB, I inverted pyramid, D inverted pyramid of C, Z inverted pyramid of Si, J inverted pyramid with APB, Fw Flat SiC(100) surface with walls, Fa Flat SiC(100) surface with aperture, SL Si(100) liquid-solid interface, FF FinFET Si (100) liquid-solid, SS Sphere Si liquid-solid" +  "\n")
        
        if Simulation == 'IN':
            file.write(path+'/'+cadfilename + "\n") # do not put comments here, only the path 
        elif Simulation == 'LA':
            file.write(path+'/'+cadfilename + "\n") # do not put comments here, only the path 
            file.write(path+'/'+tempfilename + "\n") # do not put comments here, only the path
            file.write(str(mp.Tm) + " ! Melting temperature in Kelvin" + "\n")
            file.write(str(mp.calibration_params['P0']) + " ! prefactor P0" + "\n")
            file.write(str(mp.calibration_params['A']) + " ! errf param A" + "\n")
            file.write(str(mp.calibration_params['Tflex']) + " ! errf param Tflex" + "\n")
            file.write("{:.15f}".format(LenVac) + " ! LenVac -> Top air thickness in Ang" + "\n")
            file.write("{:.15f}".format(LenNuc) + " ! LenNuc -> Nucleus radius in Ang" + "\n")
            file.write(homoflag + " ! Model Homogeneous Nucleation" + "\n")
        else:
            Sys_size=' '.join([str(s) for s in Seed_box])
            file.write(Sys_size+" ! Len1 [Len2 Len3 Len4]" + "\n")
        
        file.write(str(PtransZig)+" ! PtransZig" +  "\n")
        file.write("{:.15f}".format(mp.KMC_sf) + " ! KMC_sf -> KMC Super-Lattice parameter in Ang" + "\n")

        # Exit strategy
        file.write(ExitStrategy+" ! ExitStrategy. Iter or Time. Decides whether to use iterations or time to check simulation end and output frequency" + "\n")
        if ExitStrategy == 'Iter':
            file.write(str(IterMax)+" ! IterMax -> Max number of iterations" + "\n")
            file.write(str(OutMolMol)+" ! OutMolMol -> output frequency" + "\n")
        if ExitStrategy == 'Time':
            file.write("{:.15f}".format(TotTime) + " ! TotTime -> Max time of simulation in seconds" + "\n")
            file.write("{:.15f}".format(OutTime) + " ! OutTime -> Delta Time between consecutive print, in seconds" + "\n")

        # file.write("10000   ! exit_zeta strategy after lenz=500 " + "\n")        

        file.write(RunType+" ! RunType. R = truly random . T = not so random (use this for regression tests). C = continuation run" + "\n")
        file.write(str(IDUM)+" ! Random seed" + "\n")
        file.write(savefinalstateflag+" ! Save final state" + "\n")
        if SaveFinalState:
            file.write(path+'/'+restartfilename + "\n") # do not put comments here, only the path
        file.write(savecooflag+" ! Save final Coor" + "\n")
        if SaveCoo: # For LA this must be T
            file.write(path+'/'+coofilename + "\n") # do not put comments here, only the path

        # Probabilities
        count=0
        for indsp in range(NCrystal):
            for ind in mp.indices():
                stringa = ''.join([str(s) for s in ind])
                indd = tuple(np.insert(ind, 0, indsp))
                file.write("{} {} {} ! PtransE[{},{}] \n".format(indsp+1, stringa, PtransE[indd], indsp+1, ind))
                count+=1
        print('Number of evaporation probabilities: ', count)
        for indsp in range(NCrystal):
            for ii in range(3):
                file.write("{} {} {} ! PtransD[{},{}] \n".format(indsp+1, ii+1, PtransD[indsp,ii], indsp+1,ii+1))
        if NCov!=0:
            for indsp in range(NCrystal,NCrystal+NCov):
                for ii in range(3):
                    file.write("{} {} {} ! PtransAbs[{},{}] \n".format(indsp+1, ii+1, PtransAbs[indsp-NCrystal,ii], indsp+1,ii+1))
            for indsp in range(NCrystal,NCrystal+NCov):
                for ii in range(3):
                    file.write("{} {} {} ! PtransDes[{},{}] \n".format(indsp+1, ii+1, PtransDes[indsp-NCrystal,ii], indsp+1,ii+1))
        file.close() 
        print("DONE writing start.dat")
        return



    # ------- RUN MULSKIPS

    # Full path of the compiled mulskips executable:
    mulskipsexecutable = execpath+'/mulskips.e | tee log.txt;'
    print ("The executable directory is %s" % execpath)

    # create the directory runpath where to run the mulskips code
    # remove the runpath directory if it already exist
    try:
        shutil.rmtree(runpath)
    except OSError as e:
        print("Error: %s : %s" % (runpath, e.strerror))

    # Create new directory runpath    
    try:
        os.mkdir(runpath)
    except OSError:
        print ("Creation of the directory %s failed" % runpath)
    else:
        print ("Successfully created the directory %s " % runpath)

    # Change from current working directory path to runpath    
    try:
        os.chdir(runpath)
        print('Directory changed to: {}'.format(runpath))
    except OSError:
        print("Can't change the Current Working Directory")


    # Wrote out start.dat in runpath
    print_file_input()

    ## Save calibration parameters to txt file                       
    with open('calibration_params.txt', 'w') as ff:
        ff.write('Calibration parameters:\n')
        for key, value in mp.calibration_params.items():
            ff.write(' {} \t--> {}\n'.format(key, value))

    endok = True
    if not setup_only:
        # execute the mulskips.e program with the input file start.dat
        print('\n********************************** Starting MulSKIPs **********************************')
        make_process = subprocess.Popen(mulskipsexecutable, shell=True, stdout=subprocess.PIPE)
        while True:
            line = make_process.stdout.readline()
            if not line:break
            print('KMC:\t'+str(line.rstrip(), "utf-8")) #output to console in time
            sys.stdout.flush()
        print('********************************** Finished MulSKIPs **********************************\n')

        # !rm *.src
        # !mkdir xyz
        # !mkdir d_xyz
        # !mkdir v_xyz
        # !mkdir w_xyz
        # !mv I*_d.xyz d_xyz/
        # !mv I*_v.xyz v_xyz/
        # !mv I*_w.xyz w_xyz/
        # !mv I*.xyz xyz/

        """ Check if mulskips ended with few particles error """
        string_to_search = "few MC particles: stop MC"
        with open('log.txt', 'r') as f:
            for line in f:
                if string_to_search in line:
                    endok = False

        print('DONE MulSKIPs run. ETA: {} sec'.format(time.time()-t0))
    
    # Change from runpath to path    
    try:
        os.chdir(path)
        print('Directory changed to: {}'.format(path))
    except OSError:
        print("Can't change the Current Working Directory") 
        exit(0)


    return endok


