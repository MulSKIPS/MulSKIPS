import time
import math
import numpy as np
import os,shutil,subprocess,sys

##### --------------------------------------------------------------------------
# SETUP AND RUN MULSKIPS
##### --------------------------------------------------------------------------

###### INTRODUCI MODALITà 'INITIALIZE' E 'RUN' IN run_mulskips
### INITIALIZE: LAUNCH SOCKET E LEGGI GEOM FROM MAIN.F
### RUN: CREA START.DAt IN CURRENT FOLDER, RUN MULSKIPS, FLUSH STUFF IN RUNDIR, REMOVE STUFF FROM CURRENT FOLDER

# Function to Generate input for mulskips and run mulskips simulation 
def run_mulskips(execpath=None, runpath=None, Simulation=None, mp=None,
    PtransZig=1.0, ExitStrategy='Iter', OutMolMol=None, IterMax=None, TotTime=None, OutTime=None,
    Seed_box=[120,120,120], RunType='R', IDUM=9117116, setup_only=False,
    cadfilename=None, tempfilename=None, LenVac=None, LenNuc=None, initialState='homogeneous', LenSiGe=None,
    SaveCoo=False, coofilename=None,
    SaveFinalState=False, restartfilename=None,
    useProbTable=False, ProbTable=None, driver=None):
    """
    Generate input for mulskips and run mulskips simulation.
    """
    t0 = time.time()

    # Fundamental check
    necessary_input = [execpath, runpath, setup_only]
    if not all(v is not None for v in necessary_input):
        print('\nERROR. Please define all the necessary input (see list below):')
        print('[execpath, runpath, setup_only]')
        print(necessary_input)
        sys.exit()    



    # Fundamental check
    necessary_input = [Simulation, mp]
    if not all(v is not None for v in necessary_input):
        print('\nERROR. Please define all the necessary input (see list below):')
        print('[Simulation, mp]')
        print(necessary_input)
        sys.exit()    

    print('\n********************************** Generating MulSKIPs input file (start.dat) **********************************')
    print('*********** Input arguments:')
    print('Executable directory: ', execpath)
    print('Simulation directory: ', runpath)
    print('Type of simulation: ', Simulation)
    print('Mulskips class: ', mp)

    # Skip if MulSKIPS needs to be run in the current directory 
    path = os.getenv("PWD")
#    print(path)
#    print(os.getcwd())
#    print(os.path.abspath(runpath))

    if os.path.abspath(runpath) == path or os.path.abspath(runpath) == os.getcwd():
        print ("Running MulSKIPS in the current directory: {}".format(path))
    else:

        if runpath in ['.', './']:
            print("\nYou want to run MulSKIPS in the current directory?")
            print("Something is wrong...please set the run directory to either:\n{}\n  or this:\n{}\n".format(path, os.getcwd()))

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
    elif Simulation in ['ST']:
        necessary_input = [useProbTable]
        if not all(v is not None for v in necessary_input):
            print('\nERROR. Please define all the necessary input (see list below):')
            print('[useProbTable]')
            print(necessary_input)
            sys.exit()
        if useProbTable and ProbTable is None:
            print('\nERROR. Please define ProbTable to use in input:')
            sys.exit()
        print('Size of the system: ', Seed_box)
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
        if driver == None:
            if coofilename is None:
                print('ERROR: Please provide coofilename... ')
                sys.exit()
        else:
            print('Phases will not be written to file because sockets are ON')  

    else:
        savecooflag = 'F'

    if RunType == 'C':
        if restartfilename is None:
            print('ERROR: Please provide restartfilename... ')
            sys.exit()

    # ----- GENERATE INPUT FILE

    if Simulation != 'ST':
        if mp.__class__.__name__ == 'PhaseChange':
            
            # If Laser annealing
            if Simulation == 'LA':
                """
                Generate input for mulskips and run simulation with T map and geometry from fem.
                NB: We'll set only energetics in start.dat. Probabilities will be calculated 
                directly inside mulskips using T map from fem.
                """
                if driver==None and tempfilename==None: 
                    print('ERROR: please set \'tempfilename\' or \'driver\' flag')
                    sys.exit()
                elif driver!=None and tempfilename!=None: 
                    print('ERROR: please set \'tempfilename\' or \'driver\' flag')
                    sys.exit()

                # Check input is OK
                necessary_input = [TotTime, LenVac, LenNuc, SaveCoo, initialState, LenSiGe]
                if not all(v is not None for v in necessary_input):
                    print('\nERROR. Please define all the necessary input (see list below):')
                    print('[TotTime, LenVac, LenNuc, SaveCoo, initialState, LenSiGe]')
                    print(necessary_input)
                    sys.exit()

                if driver == None:
                    necessary_input = [tempfilename, cadfilename, coofilename, restartfilename]
                    if not all(v is not None for v in necessary_input):
                        print('\nERROR. Please define all the necessary input (see list below):')
                        print('[tempfilename, cadfilename, coofilename, restartfilename]')
                        print(necessary_input)
                        sys.exit()
                    print('Input CAD geometry file: ', path+'/'+cadfilename)    
                    print('Input temperature map file: ', path+'/'+tempfilename)  
                else:  
                    print('Input CAD and temperature files are not provided because sockets are ON')  

                print('Top air thickness : {} Ang'.format(LenVac))
                print('Nucleus radius : {} Ang'.format(LenNuc))
                print('SiGe thickness : {} Ang'.format(LenSiGe))
                
                if driver == None:
                    if SaveCoo:
                        print('Output Coor file: ', path+'/'+coofilename)    
                    else:
                        print('ERROR: Please set SaveCoo to \'T\' and provide coofilename unless sockets are used... ')
                        sys.exit()
                    print('Output checkpoint file for restart run: ', path+'/'+restartfilename)    
                
                if initialState=='homogeneous':
                    homoflag = 'T'
                elif initialState=='inhomogeneous':
                    homoflag = 'F'
                elif initialState=='tripleSF':
                    homoflag = 'D'
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
            if driver == None:
                file.write(path+'/'+cadfilename + "\n") # do not put comments here, only the path 
                file.write(path+'/'+tempfilename + "\n") # do not put comments here, only the path
            else:
                file.write('None' + "\n") # do not put comments here, only the path 
                file.write('None' + "\n") # do not put comments here, only the path 
            file.write("{} ! Melting temperature in Kelvin \n".format(' '.join([str(x) for x in mp.Tm])))
            file.write("{} ! prefactor P0 \n".format(' '.join([str(x) for x in mp.calibration_params['P0']])))
            file.write("{} ! errf param A \n".format(' '.join([str(x) for x in mp.calibration_params['A']])))
            file.write("{} ! errf param Tflex \n".format(' '.join([str(x) for x in mp.calibration_params['Tflex']])))
            file.write("{:.15f}".format(LenVac) + " ! LenVac -> Top air thickness in Ang" + "\n")
            file.write("{:.15f}".format(LenNuc) + " ! LenNuc -> Nucleus radius in Ang" + "\n")
            file.write(homoflag + " ! Initial Nucleation State [homogeneous, inhomogeneous, tripleSF] " + "\n")
            if initialState=='tripleSF':
                Sys_size=' '.join([str(s) for s in Seed_box])
                file.write(Sys_size+" ! Len1 Len2 Len3 Len4 Len5 for configuring a triple SF" + "\n")
        else:
            Sys_size=' '.join([str(s) for s in Seed_box])
            file.write(Sys_size+" ! Len1 [Len2 Len3 Len4 Len5]" + "\n")
        
        if mp.listcryZ == [14, 32]: # SiGe 
            file.write(str(mp.calibration_params['X0'][1]) + " ! Ge fraction in substrate " + "\n")
            if Simulation == 'LA':
                file.write("{:.15f}".format(LenSiGe) + " ! LenSiGe -> SiGe thickness in Ang" + "\n")
            #     file.write(str(mp.calibration_params['X'][1]) + " ! Ge fraction in liquid " + "\n")

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
            if driver == None: 
                file.write(path+'/'+coofilename + "\n") # do not put comments here, only the path
            else:
                file.write('None' + "\n") # do not put comments here, only the path

        # Probabilities
        if Simulation == 'ST' and useProbTable:
            file.write(ProbTable + "\n") # do not put comments here, only the path
        else:
            mulskipsProbThreshold = 5e-12
            count=0
            for indsp in range(NCrystal):
                for ind in mp.indices():
                    stringa = ''.join([str(s) for s in ind])
                    indd = tuple(np.insert(ind, 0, indsp))
                    if PtransE[indd] < mulskipsProbThreshold: print(f'\nWARNING: PtransE[{indd}] < mulskipsProbThreshold={mulskipsProbThreshold} --> Setting it to mulskipsProbThreshold... \n')
                    prob = np.amax([mulskipsProbThreshold, PtransE[indd]]) # ensures all probs are above minimum threshold in fortran
                    file.write("{} {} {} ! PtransE[{},{}] \n".format(indsp+1, stringa, prob, indsp+1, ind))
                    count+=1
            print('Number of evaporation probabilities: ', count)
            for indsp in range(NCrystal):
                for ii in range(3):
                    if PtransD[indsp,ii] < mulskipsProbThreshold: print(f'\nWARNING: PtransD[{indsp},{ii}] < mulskipsProbThreshold={mulskipsProbThreshold} --> Setting it to mulskipsProbThreshold... \n')
                    prob = np.amax([mulskipsProbThreshold, PtransD[indsp,ii]]) # ensures all probs are above minimum threshold in fortran
                    file.write("{} {} {} ! PtransD[{},{}] \n".format(indsp+1, ii+1, prob, indsp+1,ii+1))
            if NCov!=0:
                for indsp in range(NCrystal,NCrystal+NCov):
                    for ii in range(3):
                        if PtransAbs[indsp-NCrystal,ii] < mulskipsProbThreshold: print(f'\nWARNING: PtransAbs[{indsp-NCrystal},{ii}] < mulskipsProbThreshold={mulskipsProbThreshold} --> Setting it to mulskipsProbThreshold... \n')
                        prob = np.amax([mulskipsProbThreshold, PtransAbs[indsp-NCrystal,ii]]) # ensures all probs are above minimum threshold in fortran
                        file.write("{} {} {} ! PtransAbs[{},{}] \n".format(indsp+1, ii+1, prob, indsp+1,ii+1))
                for indsp in range(NCrystal,NCrystal+NCov):
                    for ii in range(3):
                        if PtransDes[indsp-NCrystal,ii] < mulskipsProbThreshold: print(f'\nWARNING: PtransDes[{indsp-NCrystal},{ii}] < mulskipsProbThreshold={mulskipsProbThreshold} --> Setting it to mulskipsProbThreshold... \n')
                        prob = np.amax([mulskipsProbThreshold, PtransDes[indsp-NCrystal,ii]]) # ensures all probs are above minimum threshold in fortran
                        file.write("{} {} {} ! PtransDes[{},{}] \n".format(indsp+1, ii+1, prob, indsp+1,ii+1))
        file.close() 
        print("DONE writing start.dat")
        return

    # Wrote out start.dat in runpath
    print_file_input()

    ## Save calibration parameters to txt file                       
    with open('calibration_params.txt', 'w') as ff:
        ff.write('Calibration parameters:\n')
        for key, value in mp.calibration_params.items():
            ff.write(' {} \t--> {}\n'.format(key, value))


    print('DONE MulSKIPs setup. ETA: {} sec'.format(time.time()-t0))




    # ------- RUN MULSKIPS

    # Continue only if setup_only is False
    if not setup_only:

        t0 = time.time()

        # Full path of the compiled mulskips executable:
        mulskipsexecutable = execpath+'/mulskips.e | tee log.txt;'
        print ("The executable directory is %s" % execpath)

        endok = True
        # execute the mulskips.e program with the input file start.dat
        print('\n********************************** Starting MulSKIPs **********************************')
        make_process = subprocess.Popen(mulskipsexecutable, shell=True, stdout=subprocess.PIPE)
        while True:
            line = make_process.stdout.readline()
            if not line:break
            print('KMC:\t'+str(line.rstrip(), "utf-8")) #output to console in time
            sys.stdout.flush()
        print('********************************** Finished MulSKIPs **********************************\n')


        """ Check if mulskips ended with few particles error """
        string_to_search = "few MC particles: stop MC"
        with open('log.txt', 'r') as f:
            for line in f:
                if string_to_search in line:
                    endok = False

        print('DONE MulSKIPs run. ETA: {} sec'.format(time.time()-t0))

        return endok



    # Change from runpath to path    
    try:
        os.chdir(path)
        print('Directory changed to: {}'.format(path))
    except OSError:
        print("Can't change the Current Working Directory") 
        exit(0)
    
    return


def compile_mulskips(execpath):
    # Change from runpath to path    
    try:
        os.chdir(execpath)
        print('Directory changed to: {}'.format(execpath))
    except OSError:
        print("Can't change the Current Working Directory") 
        exit(0)

    # compile the mulskips.e program 
    print('\n********************************** Starting MulSKIPS compilation **********************************')
    make_process = subprocess.Popen('make clean ; make | tee make.log;', shell=True, stdout=subprocess.PIPE)
    while True:
        line = make_process.stdout.readline()
        if not line:break
        print('MulSKIPS:\t'+str(line.rstrip(), "utf-8")) #output to console in time
        sys.stdout.flush()
    print('********************************** Finished MulSKIPS compilation **********************************\n')

    """ Check if mulskips compiled without errors """
    str1 = "Compilation terminated"
    str2 = "make: Nothing to be done for 'what'."
    endok = False
    with open('make.log', 'r') as f:
        for line in f:
            if str1 in line or str2 in line:
                endok = True

    # Change from runpath to path    
    path = os.getenv("PWD")
    try:
        os.chdir(path)
        print('Directory changed to: {}'.format(path))
    except OSError:
        print("Can't change the Current Working Directory") 
        exit(0)

    return endok


def setup_mulskips_src(execpath, lenx, leny, lenz):
    
    from tempfile import NamedTemporaryFile
    defsystemfile = execpath+'/modules/defsystem.f'
    pattern = "INTEGER, PARAMETER :: LenX="
    newstring = "       INTEGER, PARAMETER :: LenX={}, LenY={}, LenZ={}\n".format(lenx, leny, lenz)
    boxupdate, boxok = False, False
    with open(defsystemfile) as fin, NamedTemporaryFile(dir=execpath+'/modules', delete=False) as fout:
        for line in fin:
            if pattern in line:
                oldstring = line
                if line.startswith(newstring): 
                    boxok = True
                else:
                    line = newstring
                    boxupdate = True
            fout.write(line.encode('utf8'))
        os.rename(fout.name, defsystemfile)
    # Delete tmp file...
    
    # if recompile:
    #     print('Re-compiling MulSKIPS due to update of box size.')
    #     print('Old box size: {} '.format(oldstring[29:-1]))
    #     print('New box size: LenX={}, LenY={}, LenZ={}'.format(lenx, leny, lenz))
    #     endok = compile_mulskips(execpath=execpath)
    #     if not endok: 
    #         print('ERROR: Something went wrong during the compilation of MulSKIPS.')
    #         print('Further details in {}/make.log'.format(execpath))
    # elif boxok:
    #     print('Box is ok in MulSKIPS src. No need for recompiling.')
    # else:
    #     print('ERROR: Pattern not found in '+execpath+'/modules/defsystem.f')
    #     sys.exit()

    if boxupdate or boxok:
        if boxupdate:
            print('Compiling MulSKIPS with possible updated of box size.')
            print('Old box size: {} '.format(oldstring[29:-1]))
            print('New box size: LenX={}, LenY={}, LenZ={}'.format(lenx, leny, lenz))
        else:
            print('Box in MulSKIPS is ok.')
        endok = compile_mulskips(execpath=execpath)
        if not endok: 
            print('ERROR: Something went wrong during the compilation of MulSKIPS.')
            print('Further details in {}/make.log'.format(execpath))
    else:
        print('ERROR: Pattern not found in '+execpath+'/modules/defsystem.f')
        sys.exit()

    return
