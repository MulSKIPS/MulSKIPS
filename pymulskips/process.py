import time
import math
import numpy as np
import os,shutil,subprocess,sys

# ----------------------------------------------------------------------------------------------------------------------
# CLASS FOR PVD PROCESSES
# ----------------------------------------------------------------------------------------------------------------------

class PVD():

    def __init__(self, 
        substrate='SiC-3C', precursors=['Si', 'Si2C', 'SiC2'], 
        calibration_type='avrov', Tsource=1500, Tseed_center=1560,
        calibration_params=None):
        
        self.substrate=substrate
        self.precursors=precursors
        self.calibration_type=calibration_type

        # Check temperatures of source and seed (needed for computing probabilities)
        necessary_input = [Tsource, Tseed_center]
        if not all(v is not None for v in necessary_input):
            print('\nERROR. Please define all the necessary input (see list below):')
            print('[Tsource, Tseed_center]')
            print(necessary_input)
            sys.exit()
        if Tseed_center > Tsource:
            print('\nWARNING: you have chosen Tseed_center > Tsource. Proceed at your own risk...')

        self.Tsource = Tsource # K
        self.Tseed_center = Tseed_center # K

        # Some constants related to the gas phase
        self.Rgas = 8.314462618       # Gas constant J/mol/K
        self.NA = 6.02214076e23       # Avogadro Number [mol^-1]
        self.kb_ev = 8.617333262145e-5     # [eV/K]

        # By default here there is no coverage
        self.listcov = []
        self.listcovZ = []


        # Molar masses of elements (Si, C, H, Cl, B, P, Ge...)
        self.M_Si = 0.0280855  #Kg/mol
        self.M_H = 0.00100784 #Kg/mol
        self.M_C = 0.0120107 #Kg/mol
        self.M_Cl = 0.035453 #Kg/mol
        self.M_Ge = 0.07264 #Kg/mol
        self.M_B = 0.010811 #Kg/mol
        self.M_P = 0.030973762 #Kg/mol

        # Some constants related to the solid growing substrate
        if self.substrate == 'SiC-3C': 
            self.mass = self.M_Si + self.M_C # [kg/mol]
            self.rho = 3210 # kg/m^3    
            self.a0 = pow(4 * self.mass / self.NA / self.rho,1.0 / 3.0)  # [m]  #side of the cubic cell of 3C-SiC containing four Si-C dimers
            self.KMC_sf = 4.63/12.0   # KMC Super-Lattice parameter (angstrom)
            self.listcry = ['Si', 'C']
            self.listcryZ = [14, 6]
            self.listcrymass = {'Si':self.M_Si, 'C':self.M_C}

        elif self.substrate == 'Si':                                                        
            self.mass = self.M_Si # [kg/mol]
            self.rho = 2329  # [kg/m^3]
            self.a0 = pow(8 * self.mass / self.NA / self.rho,1.0 / 3.0)  # [m]  #lattice constant of Si (diamond structure, containing eight atoms per cell)
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)
            self.listcry = ['Si']
            self.listcryZ = [14]                                                      
            self.listcrymass = {'Si': self.M_Si}    
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
        self.rho_surf = self.rho * self.a0 # Kg/m^2

        # Info for the specific Process
        if self.substrate == 'SiC-3C':
            if self.precursors == ['Si', 'Si2C', 'SiC2']:
                # Some constants related to the gas-phase molecules
                M_Si = self.M_Si # [kg/mol]
                M_Si2C = 2*self.M_Si + self.M_C # [kg/mol]
                M_SiC2 = self.M_Si + 2*self.M_C # [kg/mol]
                self.precursor_masses = {'Si': M_Si, 'Si2C': M_Si2C, 'SiC2': M_SiC2}
                # Check calibration type
                if self.calibration_type is None:
                    print('ERROR: Please provide calibration_type, \'avrov\' or \'lilov\' ...')
                    sys.exit()
                # Fixed calibration for this process
                if calibration_params is None:
                    calibration_params = {
                        'fact_D_Si' : [1.10546, 0.276365, 44.2184],
                        'fact_D_C' : [200.0, 2.5, 400.0],
                        'scalef' : 0.43 # this can be the param flag to vary to calibrate
                        }
            else:
                print('ERROR: PVD process for substrate {} with precursors {} is not implemented'.format(self.substrate, self.precursors))
                sys.exit()

        elif self.substrate == 'Si':
            print('ERROR: PVD process for substrate {} is not implemented'.format(self.substrate))
            sys.exit()
    
        else:
            print('ERROR: PVD process for substrate {} is not implemented'.format(self.substrate))
            sys.exit()
    
        # Check
        calibration_params_keys = ['fact_D_Si', 'fact_D_C', 'scalef']
        for key, value in calibration_params.items():
            if key not in calibration_params_keys:
                print('  ERROR: please provide calibration parameter {} '.format(key))
                sys.exit()
        self.calibration_params = calibration_params

        print('*** Initializing PVD process for {} substrate with precursors {}'.format(self.substrate, self.precursors))
        print('precursor_masses \t--> Precursor Masses [kg/mol]:', self.precursor_masses)
        print('mass \t\t\t--> Substrate Mass [Kg/mol]:', self.mass)
        print('rho \t\t\t--> Substrate Density [Kg/m^3]:', self.rho)
        print('rho_surf \t\t--> Substrate Surface density [Kg/m^2]:', self.rho_surf)
        print('calibration_type \t--> Calibration of fluxes taken from: ', self.calibration_type)
        print('Tsource \t--> Temperature at source: ', self.Tsource)
        print('Tseed_center \t--> Temperature at seed center:', self.Tseed_center)
        print('KMC_sf \t\t\t--> Substrate KMC super-lattice constant [angstrom]:', self.KMC_sf)
        print('listcry \t\t--> Crystalline species in the substrate:', self.listcry)
        print('listcryZ \t\t--> Atomic number of Crystalline species in the substrate:', self.listcryZ)
        print('calibration_params \t--> Calibration parameters:')
        for key, value in self.calibration_params.items():
            print('     {} \t--> {}'.format(key, value))


    # -------- Some fundamental class functions to get vapor pressures and gas fluxes.
    def get_vapor_pressure(self, calibration_type,mol,temp):
        # Vapor pressures following an Arrhenius-like behaviour.
        # These are the related prefactors and exponential constants
        # for the cubic SiC from the Avrov and Lilov works.
        if calibration_type == 'avrov':
            # Avrov et al. Journal of Crystal Growth 198/199 1011-1014 (1999)
            # Avrov Arrhenius-like vapor pressures.
            A={'Si': -2.04e4, 'Si2C': -2.81e4, 'SiC2': -4.18e4}
            B={'Si': 10.82,  'Si2C': 13.28,  'SiC2': 18.18}
            Pexp=math.exp(A[mol]/temp+B[mol]) # [N/m^2]
        elif calibration_type == 'lilov':
            # Lilov et al. Cryst. Res. Technol. 28(4) 503-510 (1993)
            # Lilov  1 1500-2000 Arrhenius-like vapor pressures.
            A_1={'Si':-27499.8,'Si2C':-34177.2,'SiC2':-34075.8}      # Temperature range 1500-2000 [K]
            A_2={'Si':-27261.32,'Si2C':-33500.33,'SiC2':-33526.61}   # Temperature range 2000-2546 [K]
            B_1={'Si':12.8114,'Si2C':15.1877,'SiC2':15.4274}         # Temperature range 1500-2000 [K]
            B_2={'Si':12.6921,'Si2C':14.8493,'SiC2':15.1528}         # Temperature range 2000-2546 [K]
            if 1500 <= temp <=2000:
                Pexp=math.exp(A_1[mol]/temp+B_1[mol]) # [N/m^2]
            # Lilov  2 2000-2546 Arrhenius-like vapor pressures.
            elif 2000 < temp <=2546:
                Pexp=math.exp(A_2[mol]/temp+B_2[mol]) # [N/m^2]
            else:
                print("error, temperature range not yet implemented.")
        else:
            print("error, calibration_type wrong.")
        return Pexp     

    def get_flux(self, calibration_type, mol, temp):
        p_dep = self.get_vapor_pressure(calibration_type,mol,temp)   # Pascal
        coeff = pow(2.*math.pi*self.precursor_masses[mol]*self.Rgas*temp,-0.5)
        flux = p_dep * coeff # mol/s/m^2 
        return flux

    def get_fluxes(self, calibration_type, temp):
        # This function give the <data> dictionary of "Supersaturation SiC$_2$", "Supersaturation Si$_2$C"
        # and "3C-SiC Growth rate [micron/hour]" for the Avrov and Lilov calibration of the vapor pressuare within
        # the growth chamber.
        fluxes = {}
        for mol in self.precursors:
            fluxes[mol] = self.get_flux(calibration_type,mol,temp) # mol/s/m^2
        return fluxes

    def get_fluxes_evap(self, loc_fluxes):
        # This function gives the fluxes for evaporation 
        # To Do: Define stoichiometry vector (e.g. v = [1.0, 2.0, 1.0]) above in init and then 
        # here do a scalar product e.g. 
        # fv = np.fromiter(loc_fluxes.values(), dtype=float)
        # Si_flux_ev =  np.dot(v, fv))
        # In this way I can avoid the if checks below...
        # Also, this can be a simple way of including reaction energetics. 
        # For example, if Si2C does not contribute to Si_flux_ev with 2 Si atoms one 
        # could use the vector [1.0, 1.8, 1.0] as v .
        if self.substrate == 'SiC-3C':
            if self.precursors == ['Si', 'Si2C', 'SiC2']:
                Si_flux_ev = loc_fluxes['Si'] + 2.0*loc_fluxes['Si2C'] + loc_fluxes['SiC2']
                C_flux_ev = loc_fluxes['Si2C'] + 2.0*loc_fluxes['SiC2']
                return Si_flux_ev, C_flux_ev
            else:
                print('ERROR: not implemented')
        else:
            print('ERROR: not implemented')

    # -------------- AUXILIARY ROUTINES 
    # Indices is not used inside the class, but it's useful to write out start.dat (can be fixed better...)
    def indices(self):
        """
        Returns a numpy array containing the indices of elements in a probability matrix for a system with:  
        ncry: number of crystalline species 
        ncov: number of coverage species
        """
        ncry, ncov = len(self.listcry), len(self.listcov)
        x=[0,1,2,3]
        # Find all the possible permutations of elements in x
        import itertools as it
        perm = np.asarray(list(it.product(x, repeat=ncry+ncov)))
        # Exclude those where the sum of all elements is larger than 4 (i.e. the total number of nearest neighbours)
        sumperm = np.sum(perm, axis=1)
        todelete = np.where(sumperm>4)[0]
        # Exclude [0,0,0,0]
        todelete = np.append(todelete, np.where(sumperm==0)[0])
        # Exclude those where the first ncry elements are zero (there is always at least 1 crystal neighbour)
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)==0)[0])
        # Exclude those where the sum of the first ncry elements is larger than 3 (a crystal site cannot have more than 3 crystal neighbours in a tetrahedral configuration) 
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)>3)[0])
        # Array of survivors
        inds = np.delete(perm, todelete, axis=0)
        return inds

    def get_dims(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        return np.ones(ncry+ncov, dtype=int)*4    


    # -------------- SETUP EVAPORATION FREQUENCIES
    def get_energetics_evap(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        Etot = np.zeros((ncry,*dims))  # energetics of evaporation in eV
        
        if self.substrate == 'SiC-3C':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int)*4
            E0 = np.zeros((ncry,*dims0))

            f=4 # factor for the energetics interpolation
            # Coordination 1 of Si from 111 C-terminated 4x1x2 bigdft surface
            E0[0,0,1] = 6.83 # -4.39 BigDFT # -6.83 # QE 
            # Coordination 2 of Si from 100 Si-terminated 3x1x3 dimer surface
            E0[0,0,2] = 7.04
            # Coordination 3 of Si from 111 Si-terminated 4x1x2 surface
            E0[0,0,3] = 7.77
            # Coordination 1 of Si adatom optimized on 111 Si-terminated 4x1x2 surface
            E0[0,1,0] = 6.72
            E0[0,1,1] = (f*E0[0,1,0]+E0[0,0,2])/(f+1) # -6.8 # E0[0,1,0] < E0[0,1,1] < E0[0,0,2]
            E0[0,1,2] = (f*E0[0,0,2]+E0[0,0,3])/(f+1) # -6.8 # E0[0,0,2] < E0[0,1,2] < E0[0,0,3]
            E0[0,2,1] = (f*E0[0,1,1]+E0[0,0,3])/(f+1) # -6.8 # E0[0,1,1] < E0[0,2,1] < E0[0,0,3]
            E0[0,2,0] = (f*E0[0,1,0]+E0[0,2,1])/(f+1) # -6.8 # E0[0,1,0] < E0[0,2,0] < E0[0,2,1]
            E0[0,3,0] = (f*E0[0,2,0]+E0[0,2,1])/(f+1) # -6.8 # E0[0,2,0] < E0[0,3,0] < E0[0,2,1]
            # Coordination 1 of C from 111 Si-terminated 4x1x2 surface
            E0[1,1,0] = 9.10
            # Coordination 2 of C from 100 C-terminated 3x1x3 dimer surface
            E0[1,2,0] = 11.43
            # Coordination 3 of C from 111 C-terminated 4x1x2 bigdft surface
            E0[1,3,0] = 11.80
            # Coordination 1 of C adatom optimized on 111 C-terminated 4x1x2 bigdft surface
            E0[1,0,1] = 7.76
            E0[1,1,1] = (f*E0[1,1,0]+E0[1,2,0])/(f+1) #  -9.8 # E0[1,1,0] < E0[1,1,1] < E0[1,2,0]
            E0[1,1,2] = (f*E0[1,1,1]+E0[1,3,0])/(f+1) # -10.0 # E0[1,1,1] < E0[1,1,2] < E0[1,3,0]
            E0[1,0,2] = (f*E0[1,0,1]+E0[1,1,2])/(f+1) # -8.1 # E0[1,0,1] < E0[1,0,2] < E0[1,1,2] 
            E0[1,2,1] = (f*E0[1,2,0]+E0[1,3,0])/(f+1) # -11.6 # E0[1,2,0] < E0[1,2,1] < E0[1,3,0]
            E0[1,0,3] = (f*E0[1,0,2]+E0[1,1,2])/(f+1) # -9.1 # E0[1,0,2] < E0[1,0,3] < E0[1,1,2]
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()

        for i in range(ncry):
            Etot[i] = E0[i]

        return Etot

    def get_PtransE(self):    

        # Set relevant temperature
        Temperature = self.Tseed_center

        # Energetics for evaporation
        ene = self.get_energetics_evap()
        print('\nEvaporation energetics:\n {}'.format(ene))

        # Prefactor
        freq_0 = self.mass / self.rho_surf # m^2 / mol
        # freq_0 = self.NA / self.at_density  # m^2/mol

        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        PtransEtot = np.zeros((ncry,*dims))
        Bolt = np.zeros((ncry,*dims))                

        loc_fluxes = self.get_fluxes(self.calibration_type, Temperature) # mol/s/m^2
        fluxes_evap = self.get_fluxes_evap(loc_fluxes)

        for i in range(ncry):
            ll = np.zeros(ncry, dtype=int)
            ll[-(1+i)]=2
            # if ncry==1:
            #     ll[-(1+i)]=1    
            # elif ncry>1:
            #     ll[-(1+i)]=2    
            enecoo2 = ene[tuple(np.r_[i,ll])]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            Bolt[i]=np.exp(-(ene[i] - enecoo2)/(self.kb_ev * Temperature))
            PtransEtot[i] = freq_0 * fluxes_evap[i] * Bolt[i] # at/s

        return PtransEtot

    # -------------- SETUP DEPOSITION FREQUENCIES
    def get_PtransD(self):    

        # Set relevant temperature
        Temperature = self.Tsource

        # Energetics for evaporation
        ncry = len(self.listcry)
        ene = np.zeros((ncry,3))
        ene_evap = self.get_energetics_evap()
        for i in range(ncry):
            ll = np.zeros(ncry, dtype=int)
            ll[-(1+i)]=2
            # if ncry==1:
            #     ll[-(1+i)]=1    
            # elif ncry>1:
            #     ll[-(1+i)]=2    
            ene[i] = ene_evap[tuple(np.r_[i,ll])]
        print('\nDeposition energetics:\n {}'.format(ene))

        # Prefactor
        freq_0 = self.mass / self.rho_surf # m^2 / mol
        # freq_0 = self.NA / self.at_density  # m^2/mol

        # Calculate Frequencies
        PtransDtot = np.zeros((ncry,3))
        Bolt = np.zeros((ncry,3))

        loc_fluxes = self.get_fluxes(self.calibration_type, Temperature) # mol/s/m^2
        fluxes_evap = self.get_fluxes_evap(loc_fluxes)
        
        for i in range(ncry):
            ll = np.zeros(ncry, dtype=int)
            ll[-(1+i)]=2
            # if ncry==1:
            #     ll[-(1+i)]=1    
            # elif ncry>1:
            #     ll[-(1+i)]=2
            enecoo2 = ene[i,1]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            Bolt[i]=np.exp(-(ene[i] - enecoo2)/(self.kb_ev * Temperature))
            PtransDtot[i] = freq_0 * fluxes_evap[i] * Bolt[i] \
            * self.calibration_params['scalef'] \
            * self.calibration_params[f'fact_D_{self.listcry[i]}'] # at/s

        return PtransDtot



# ----------------------------------------------------------------------------------------------------------------------
# CLASS FOR CVD PROCESSES
# ----------------------------------------------------------------------------------------------------------------------

class CVD():

    def __init__(self, expdata, calibration_params=None, miller=100):  
        
        print('Reading input from provided dictionary')
        self.expdata = expdata

        self.substrate = expdata['substrate'] 
        self.precursors = expdata['precursors'] 

        self.temp = expdata['temperature']
        self.pressure = {}
        for mol in self.precursors:
            self.pressure[mol] = expdata['pressure_{}'.format(mol)]

        self.miller = miller

        # Some constants related to the gas phase
        self.Rgas = 8.314462618       # Gas constant J/mol/K
        self.NA = 6.02214076e23       # Avogadro Number [mol^-1]
        self.kb_ev = 8.617333262145e-5     # [eV/K]
        self.h = 4.135667696e-15 # eV*s

        # Molar masses of elements (Si, C, H, Cl, B, P, Ge...)
        self.M_Si = 0.0280855  #Kg/mol
        self.M_H = 0.00100784 #Kg/mol
        self.M_C = 0.0120107 #Kg/mol
        self.M_Cl = 0.035453 #Kg/mol
        self.M_Ge = 0.07264 #Kg/mol
        self.M_B = 0.010811 #Kg/mol
        self.M_P = 0.030973762 #Kg/mol

        # Some constants related to the solid growing substrate
        if self.substrate == 'SiC-3C': 
            self.mass = self.M_Si + self.M_C # [kg/mol]
            self.rho = 3210 # kg/m^3    
            self.listcry = ['Si', 'C']
            self.listcryZ = [14, 6]
            self.a0 = pow(4*self.mass/self.NA/self.rho,1.0/3.0) # [m] side of the cubic conventional cell of 3C-SiC containing four Si-C dimers
            self.KMC_sf = 4.63/12.0   # KMC Super-Lattice parameter (angstrom)
        elif self.substrate == 'Si': 
            self.mass = self.M_Si # [kg/mol]
            self.rho = 2329  # [kg/m^3]
            self.a0 = pow(8 * self.mass / self.NA / self.rho,1.0 / 3.0)  # [m]  #lattice constant of Si (diamond structure, containing eight atoms per cell)
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)
            self.listcry = ['Si']
            self.listcryZ = [14]
            self.listcrymass = {'Si':self.M_Si}
            self.Tm = 1688 # K # Melting
            self.Td = 640 # K # Debye #corrected by 0.6 to estimate SURFACE Td
            if self.miller == 100:
                self.at_density = 6.78e18 # atoms/m^2
            elif self.miller == 110:
                self.at_density = 9.59e18 # atoms/m^2
            elif self.miller == 111:
                self.at_density = 7.83e18 # atoms/m^2
            else:
                print('ERROR: Surface with Miller indices {} is not implemented'.format(self.miller))

        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
        # Below is an approiximated rho_surf. For Miller indices differientiation, we should include them for each face in the lines above 
        self.rho_surf = self.rho * self.a0 # Kg/m^2

        # Info for the specific Process    
        if self.substrate == 'Si':
            if self.precursors == ['SiH4', 'H2']:
                # Some constants related to the gas-phase molecules
                M_SiH4 = self.M_Si + 4*self.M_H # [kg/mol]
                M_H2 = 2*self.M_H # [kg/mol]
                self.precursor_masses = {'SiH4': M_SiH4, 'H2': M_H2}
                self.listcov = ['H']
                self.listcovZ = [1]

                # Check
                if calibration_params is None:
                    print('ERROR: Please provide calibration parameters dictionary...')
                    sys.exit()
                calibration_params_keys = ['Ed1', 'Ed2', 'Ed3', 'Ee1', 'Ee2', 'Ee3', 'Eabs1', 'Eabs2', 'Eabs3', 
                    'Edes1', 'Edes2', 'Edes3', 'deltaH', 'Ereact_SiH4', 'Ereact_H2', 
                    'k_Si_SiH4', 'k_H_SiH4', 'k_H_H2', 'scalef', 
                    'aads', 'bads', 's0', 'Ades', 'Aev', 'Bev', 'Adep']
                for key, value in calibration_params.items():
                    if key not in calibration_params_keys:
                        print('  ERROR: please provide calibration parameter {} '.format(key))
                        sys.exit()

            else:
                print('ERROR: CVD process for substrate {} with precursors {} is not implemented'.format(self.substrate, self.precursors))
                sys.exit()
        else:
            print('ERROR: CVD process for substrate {} is not implemented'.format(self.substrate))
            sys.exit()

        # Assign calibration params
        self.calibration_params = calibration_params
                
        print('*** Initializing CVD process for {} substrate with precursors {}'.format(self.substrate, self.precursors))
        print('mass \t\t\t--> Substrate Mass [Kg/mol]:', self.mass)
        print('rho \t\t\t--> Substrate Density [Kg/m^3]:', self.rho)
        print('rho_surf \t\t--> Substrate Surface density [Kg/m^2]:', self.rho_surf)
        print('at_density \t\t--> Substrate Surface atomic density [at/m^2]:', self.at_density)
        print('temp \t\t\t--> Substrate temperature [K]:', self.temp)
        print('precursor_masses \t--> Precursor Masses [kg/mol]:', self.precursor_masses)
        print('pressure \t\t--> Precursors pressure [Pa]:', self.pressure)
        print('listcry \t\t--> Crystalline species in the substrate:', self.listcry)
        print('listcryZ \t\t--> Atomic number of Crystalline species in the substrate:', self.listcryZ)
        print('listcov \t\t--> Coverage species:', self.listcov)
        print('listcovZ \t\t--> Atomic number of coverage species:', self.listcovZ)
        print('KMC_sf \t\t\t--> Substrate KMC super-lattice constant [angstrom]:', self.KMC_sf)
        print('calibration_params \t--> Calibration parameters:')
        for key, value in self.calibration_params.items():
            print('     {} \t--> {}'.format(key, value))

    # Some fundamental class functions to get vapor pressures and gas fluxes.
    # # Some fundamental class functions to get vapor pressures and gas fluxes.

    # def get_vapor_pressure(self, mol): 
    #     if mol == 'Si': 
    #         # vapor pressure of Si
    #         self.pvA=23877
    #         self.pvB=12.19
    #         Pv = pow(10, -self.pvA/self.temp + self.pvB)  # [Pascal]       # PROVVISORIA: calcolata da log(PSi) [Pa] = –18558/T + 11.73 a T=1593K
    #         # print('Using this vapor pressure for si: {}'.format(Pv))
    #     elif mol == 'H': 
    #         # vapor pressure of Si
    #         Pv = 55555555555 # [Pascal]
    #     else:
    #         print('ERROR: Vapor pressure for {} is not implemented'.format(mol))
    #         sys.exit()
    #     return Pv

    def get_flux(self,mol): #,vapor_pressure=False):
        # if vapor_pressure:
        #     p_dep = self.get_vapor_pressure(mol)  # [Pascal]
        #     coeff = pow(2.*math.pi*self.mass*self.Rgas*self.temp,-0.5) 
        # else:
        #     p_dep = self.pressure[mol]   # Pascal
        #     coeff = pow(2.*math.pi*self.precursor_masses[mol]*self.Rgas*self.temp,-0.5)
        p_dep = self.pressure[mol]   # Pascal
        coeff = pow(2.*math.pi*self.precursor_masses[mol]*self.Rgas*self.temp,-0.5)
        flux = p_dep * coeff # mol/s/m^2 
        return flux

    def get_fluxes(self):
        fluxes = {}
        for mol in self.precursors:
            fluxes[mol] = self.get_flux(mol) # mol/s/m^2
        return fluxes



    # -------------- SETUP DEPOSITION FREQUENCIES

    def get_energetics_depo(self):
        ncry = len(self.listcry)
        Etot=np.zeros((ncry,3))  # energetics of deposition in eV
        
        E0=np.zeros((ncry,3)) 
        if self.substrate == 'Si':
            E0[0,0]= self.calibration_params['Ed1'] # coor 1 
            E0[0,1]= self.calibration_params['Ed2'] # coor 2
            E0[0,2]= self.calibration_params['Ed3'] # coor 3 
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
            
        # Set E0 in Etot 
        for i in range(ncry):
            Etot[i] = E0[i]
        
        return Etot

    def get_fluxes_depo(self):
        
        ncry, ncov = len(self.listcry), len(self.listcov)

        # Get precursors fluxes on surface
        loc_fluxes = self.get_fluxes() # mol/s/m^2
        
        # Get atomic fluxes on surface  
        tot_fluxes = np.zeros(ncry+ncov)
        if self.substrate == 'Si':
            # To Do: Define stoichiometry vector (e.g. v = [1.0, 2.0, 1.0]) above in init and then 
            # here do a scalar product e.g. 
            # fv = np.fromiter(loc_fluxes.values(), dtype=float)
            # Si_flux_ev =  np.dot(v, fv))
            # In this way I can avoid the if checks below...
            # Also, this can be a simple way of including reaction energetics. 
            # For example, if Si2C does not contribute to Si_flux_ev with 2 Si atoms one 
            # could use the vector [1.0, 1.8, 1.0] as v .
            if self.precursors == ['SiH4', 'H2']:
                # Consider the reaction barrier
                reactf_SiH4 = np.exp(- self.calibration_params['Ereact_SiH4'] /(self.kb_ev * self.temp)) 
                reactf_H2 = np.exp(- self.calibration_params['Ereact_H2'] /(self.kb_ev * self.temp)) 
                tot_fluxes[0] = loc_fluxes['SiH4'] * self.calibration_params['k_Si_SiH4'] * reactf_SiH4  # Silicon
                tot_fluxes[1] = loc_fluxes['SiH4'] * self.calibration_params['k_H_SiH4'] * reactf_SiH4 + \
                                loc_fluxes['H2'] * self.calibration_params['k_H_H2'] * reactf_H2  # Hydrogen
            else:
                print('ERROR: not implemented'); sys.exit()
        else:
            print('ERROR: Not implemented'); sys.exit()

        return tot_fluxes

    def get_PtransD(self):
        
        # Energetics for deposition
        ene = self.get_energetics_depo()
        print('\nDeposition energetics:\n {}'.format(ene))

        # Prefactor
        #freq_0 = self.mass / self.rho_surf # m^2 / mol
        freq_0 = self.NA / self.at_density  # m^2/mol

        # Fluxes
        fluxes_depo = self.get_fluxes_depo()

        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        PtransDtot = np.zeros((ncry,3))
        Bolt = np.zeros((ncry,3))
        for i in range(ncry):
            ll = np.zeros(ncry,dtype=int)
            if ncry==1:
                ll[-(1+i)]=1
            elif ncry>1:
                ll[-(1+i)]=2
            enecoo2 = ene[tuple(np.r_[i,ll])]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            #print('Check for Coor 2 energy used to normalize (eV): ', enecoo2)
            # Bolt[i] = np.exp(-(ene[i] - enecoo2) / (self.kb_ev * self.temp))
            # print('Check for Coor 2 energy used to normalize (eV): ', enecoo2)
            Bolt[i]=np.exp(-ene[i]/(self.kb_ev * self.temp))
            PtransDtot[i] = self.calibration_params['scalef']*self.calibration_params['Adep'] * freq_0 * fluxes_depo[i] * Bolt[i]  # at/s

        return PtransDtot


    # -------------- SETUP EVAPORATION FREQUENCIES

    # -------------- AUXILIARY ROUTINES
    def indices(self):
        """
        Returns a numpy array containing the indices of elements in a probability matrix for a system with:  
        ncry: number of crystalline species 
        ncov: number of coverage species
        """
        ncry, ncov = len(self.listcry), len(self.listcov)
        x=[0,1,2,3]
        # Find all the possible permutations of elements in x
        import itertools as it
        perm = np.asarray(list(it.product(x, repeat=ncry+ncov)))
        # Exclude those where the sum of all elements is larger than 4 (i.e. the total number of nearest neighbours)
        sumperm = np.sum(perm, axis=1)
        todelete = np.where(sumperm>4)[0]
        # Exclude [0,0,0,0]
        todelete = np.append(todelete, np.where(sumperm==0)[0])
        # Exclude those where the first ncry elements are zero (there is always at least 1 crystal neighbour)
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)==0)[0])
        # Exclude those where the sum of the first ncry elements is larger than 3 (a crystal site cannot have more than 3 crystal neighbours in a tetrahedral configuration) 
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)>3)[0])
        # Array of survivors
        inds = np.delete(perm, todelete, axis=0)
        return inds

    def deltaEcov(self, indici): # perturbation for evaporation energies
        """
        If ncov is zero, return zero
        If ncov is zero, return deltaE: energy barrier for evaporation in eV, induced by coverage
        Note:
        H --> deltaE > 0
        Cl --> deltaE < 0
        BUT this is only temporary!
        """
        # Optionally provide an extra adjustment 
        ncry, ncov = len(self.listcry), len(self.listcov)
        extra = 1  
        if self.substrate == 'SiC-3C':
            extra = 1e1
        
        deltaH = self.calibration_params['deltaH'] # eV     # CALIBRATE
        deltaCl = 1e-2 # eV    # CALIBRATE
        
        deltaE=0
        if ncov!=0:
            covn = indici[ncry:]
            for ii,nn in enumerate(covn):
                # CASE: H
                if self.listcov[ii] == 'H': 
                    deltaE += nn*deltaH*extra  
                # CASE: Cl
                elif self.listcov[ii] == 'Cl':
                    deltaE -= nn*deltaCl*extra 
                else:
                    print('ERROR: Coverage species {} is not implemented'.format(self.listcov[ii]))         
                    exit(0)
        return deltaE

    def get_dims(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        return np.ones(ncry+ncov, dtype=int)*4    


    def get_energetics_evap(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        Etot=np.zeros((ncry,*dims))  # energetics of evaporation in eV
        
        if self.substrate == 'Si':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int)*4
            E0 = np.zeros((ncry,*dims0))
            E0[0,1] = self.calibration_params['Ee1'] # 1coor, no coverage
            E0[0,2] = self.calibration_params['Ee2'] # 2coor, no coverage 
            E0[0,3] = self.calibration_params['Ee3'] # 3coor, no coverage
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
            
        # Set E0 in Etot and allocate extra dimensions for coverages
        if ncov==0:
            for i in range(ncry):
                Etot[i] = E0[i]
        elif ncov==1:
            for i in range(ncry):
                Etot[i,...,0] = E0[i]
        elif ncov==2:
            for i in range(ncry):
                Etot[i,...,0,0] = E0[i]
        else:
            print('ERROR: ncov too high, Not implemented')
            sys.exit()
                
        # Set effect of coverages in Etot
        for indsp in range(ncry): # loop over species
            for ind in self.indices(): # loop over elements in the matrix
                if sum(ind[ncry:])!=0: # if there is at least one coverage neighbour
                    indd = tuple(np.r_[indsp,ind]) # general index format in Etot
                    indd0 = tuple(np.r_[indsp,ind[:ncry],np.zeros(ncov, dtype=int)]) # indices for submatrix with 0 coverages 
                    Etot[indd] = Etot[indd0] + self.deltaEcov(ind) 
                    # the minus sign is because the Boltzmann exponential will be positive, 
                    # and we want evap probability to decrease if deltaE > 0
        
        return Etot

    def get_PtransE(self):    

        # Energetics for evaporation
        ene = self.get_energetics_evap()
        print('\nEvaporation energetics:\n {}'.format(ene))

        # Prefactor
        #freq_0 = self.mass / self.rho_surf # m^2 / mol
        freq_0 = self.NA / self.at_density  # m^2/mol
        
        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        PtransEtot = np.zeros((ncry,*dims))
        Bolt, norm = np.zeros((ncry,*dims)), np.zeros(ncry)                
        fluxes_from_vapor_pressure = np.zeros(ncry)
        for i in range(ncry):
            ll = np.zeros(ncry,dtype=int)
            ll[-(1+i)]=2
            enecoo2 = ene[tuple(np.r_[i,ll,np.zeros(ncov, dtype=int)])]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            # Bolt[i]=np.exp(-(ene[i] - enecoo2)/(self.kb_ev * self.temp))
            # print('Check for Coor 2 energy used to normalize (eV):', enecoo2)
            Bolt[i] = np.exp(-ene[i]/ (self.kb_ev * self.temp))
            # Get fluxes from vapor pressure
            # fluxes_from_vapor_pressure[i] = self.get_flux(self.listcry[i], vapor_pressure=True) # mol/s/m^2
            # PtransEtot[i] = freq_0*fluxes_from_vapor_pressure[i]*Bolt[i]* self.calibration_params['scalef'] * self.calibration_params['Aev'] # at/s
            PtransEtot[i] = self.calibration_params['scalef'] * (self.Td * self.kb_ev / self.h) * (self.calibration_params['Aev'] * self.temp + self.calibration_params['Bev']) * Bolt[i] # at/s

        return PtransEtot

    # -------------- SETUP ABSORPTION FREQUENCIES

    def get_energetics_abs(self):
        ncov = len(self.listcov)
        Etot=np.zeros((ncov,3))  # energetics of absorption in eV
        
        if self.listcov == ['H']:
            E0=np.zeros((ncov,3)) 
            E0[0,0]= self.calibration_params['Eabs1'] # coor 1   # <---- TO BE CALIBRATED
            E0[0,1]= self.calibration_params['Eabs2'] # coor 2   # <---- TO BE CALIBRATED
            E0[0,2]= self.calibration_params['Eabs3'] # coor 3   # <---- TO BE CALIBRATED
        else:
            print('ERROR: Not implemented')
            exit(0)
            
        # Set E0 in Etot 
        for i in range(ncov):
            Etot[i] = E0[i]
        
        return Etot

    def get_PtransAbs(self):    
        # Initialize probability matrix with dimension depending on NCrystal and NCov
        
        # Energetics for absorption
        ene = self.get_energetics_abs()
        print('\nAbsorption energetics:\n {}'.format(ene))

        # Prefactor
        #freq_0 = self.mass / self.rho_surf # m^2 / mol
        freq_0 = self.NA/self.at_density #m^2/mol
        
        # Fluxes
        fluxes_depo = self.get_fluxes_depo()

        # Energy barrier(T) #MMonCa
        self.eads=self.calibration_params['aads']-self.calibration_params['bads']*self.kb_ev*self.temp

        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        PtransAbstot = np.zeros((ncry,3))
        Bolt = np.zeros((ncry,3))
        for i in range(ncov):
            ll = np.zeros(ncov,dtype=int)
            if ncry==1:
                ll[-(1+i)]=1    
            elif ncry>1:
                ll[-(1+i)]=2
            # enecoo2 = ene[tuple(np.r_[i,ll])]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            #print('Check for Coor 2 energy used to normalize (eV): ', enecoo2)
            #Bolt[i]=np.exp(-(ene[i] - enecoo2)/(self.kb_ev * self.temp))
            Bolt = np.exp(-self.eads/ (self.kb_ev * self.temp))
            PtransAbstot[i] = self.calibration_params['scalef'] * freq_0 * self.calibration_params['s0'] * fluxes_depo[ncry+i] * Bolt # at/s
            #PtransAbstot[i] = fluxes_depo[ncry+i] #piecewise check
            #PtransAbstot[i] = self.pressure['H2']*pow(2*math.pi*self.M_H*2*self.Rgas*self.temp,-0.5) #piecewise check

        return PtransAbstot 

    # -------------- SETUP DESORPTION FREQUENCIES

    def get_energetics_des(self):
        ncov = len(self.listcov)
        Etot=np.zeros((ncov,3))  # energetics of absorption in eV
        
        if self.listcov == ['H']:
            E0=np.zeros((ncov,3)) 
            E0[0,0]= self.calibration_params['Edes1'] # coor 1   # <--- TO BE CALIBRATED
            E0[0,1]= self.calibration_params['Edes2'] # coor 2   # <--- TO BE CALIBRATED
            E0[0,2]= self.calibration_params['Edes3'] # coor 3   # <--- TO BE CALIBRATED
        else:
            print('ERROR: Not implemented')
            exit(0)
            
        # Set E0 in Etot 
        for i in range(ncov):
            Etot[i] = E0[i]
        
        return Etot

    def get_PtransDes(self):    
        # Initialize probability matrix with dimension depending on NCrystal and NCov

        # Energetics for absorption
        ene = self.get_energetics_des()
        print('\nDesorption energetics:\n {}'.format(ene))

        # Prefactor
        # we have no fluxes now, so we should convert freq in s^-1
        freq = self.kb_ev*self.Td / self.h   # clock of surface
        print('freq: ',freq)

        ncry, ncov = len(self.listcry), len(self.listcov)
        PtransDestot = np.zeros((ncov,3))
        Bolt = np.zeros((ncov,3))
        for i in range(ncov):
            Bolt[i]=np.exp(-ene[i]/(self.kb_ev * self.temp))   # non normalizziamo ocn enecoo2 qui, perché quel ragionamento sul coor=2 che riproduce l'andamento medio vale quando il prefattore è funzione del flusso
            PtransDestot[i] = self.calibration_params['scalef'] * self.calibration_params['Ades'] * Bolt[i] # at/s

        return PtransDestot


# ----------------------------------------------------------------------------------------------------------------------
# CLASS FOR LA PROCESSES
# ----------------------------------------------------------------------------------------------------------------------

class PhaseChange(): # when we have Ge and SiGE calibration, rename this as melting_MulskipsProcess()

    def __init__(self, material, temp=None, calibration_params=None):  
        
        self.substrate = material
        if temp is not None:
            self.temp = temp

        # Some constants related to the gas phase
        self.Rgas = 8.314462618       # Gas constant J/mol/K  # 8.314462618*6.242e18 eV/mol/K 
        self.NA = 6.02214076e23       # Avogadro Number [mol^-1]
        self.kb_ev = 8.617333262145e-5     # [eV/K]

        # By default here there is no coverage
        self.listcov = []
        self.listcovZ = []

        # Some constants related to the solid growing substrate
        if self.substrate == 'Si': 
            self.mass = 0.0280855 # [kg/mol]
            self.rho = 2329 # kg/m^3
            self.listcry = ['Si']
            self.listcryZ = [14]
            self.Tm = 1688 # K
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)

            if calibration_params is None:
                delta=0.03 #ev
                ED2=0.96 #eV #from article
                calibration_params = {
                    # Energy barrier for deposition of solid Si in a site with coor 1, 2 or 3 (PsiL-PsiS in the article...) 
                    'Ed1': ED2-delta, 'Ed2': ED2, 'Ed3': ED2+delta,  #Ed2 is fixed, vary the others
                    # Prefactor 
                    'P0': 1.33e17, #1.56e17
                    # constant in damping exp for depo probs
                    'A': 280, #settato per error function #350
                    'Tflex': 1080 # un po' meno della T del flesso 
                    }

        elif self.substrate == 'Ge':
            self.mass = 0.07263 # [kg/mol]
            self.rho = 2329 + 3493 * 1- 499 * pow(1, 2)  # kg/m^3
            self.listcry = ['Ge']
            self.listcryZ = [32]
            self.Tm = 1210  # K ! Change depending on Ge concentration
            self.KMC_sf = 5.66 / 12.0  # KMC Super-Lattice parameter (angstrom)
        elif self.substrate == 'SiGe02fake':
            self.mass = 0.07263*0.2 + 0.0280855*0.8  # [kg/mol]
            self.rho = 2329+3493*0.2-499*pow(0.2,2)  # kg/m^3
            self.listcry = ['SiGe02fake']
            self.listcryZ = [14] #[14*0.8+32*0.2] #deve essere un intero
            self.Tm = 1210*0.2 + 1688*(1.0-0.2) # K ! Change depending on Ge concentration
            self.KMC_sf = (5.43*0.8+5.66*0.2)/12.0  # KMC Super-Lattice parameter (angstrom)
        elif self.substrate == 'SiGe04fake':
            self.mass = 0.07263*0.4 + 0.0280855*0.6  # [kg/mol]
            self.rho = 2329+3493*0.4-499*pow(0.4,2)  # kg/m^3
            self.listcry = ['Si04fake']
            self.listcryZ = [14] #[14*0.6+32*0.4]
            self.Tm = 1210*0.4 + 1688*(1-0.4)  # K ! Change depending on Ge concentration
            self.KMC_sf = (5.43*0.6+5.66*0.4)/ 12.0  # KMC Super-Lattice parameter (angstrom)
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
        a0 = pow(8*self.mass/self.NA/self.rho,1.0/3.0) # [m] lattice parameter for crystalline Si (diamond structure, 8 atoms per cell)
        self.rho_surf = self.rho * a0 # Kg/m^2
                
        # Check
        calibration_params_keys = ['Ed1', 'Ed2', 'Ed3', 'P0', 'A', 'Tflex']
        for key, value in calibration_params.items():
            if key not in calibration_params_keys:
                print('  ERROR: please provide calibration parameter {} '.format(key))
                sys.exit()
        self.calibration_params = calibration_params

        print('*** Initializing melting/solidification process for {} '.format(self.substrate))
        print('mass \t\t\t--> Substrate Mass [Kg/mol]:', self.mass)
        print('rho \t\t\t--> Substrate Density [Kg/m^3]:', self.rho)
        print('rho_surf \t\t--> Substrate Surface density [Kg/m^2]:', self.rho_surf)
        if temp is not None:
            print('temp \t\t\t--> Substrate temperature [K]:', self.temp)
        print('Tm \t\t\t--> Substrate melting temperature [K]:', self.Tm)
        print('KMC_sf \t\t\t--> Substrate KMC super-lattice constant [angstrom]:', self.KMC_sf)
        print('listcry \t\t--> Crystalline species in the substrate:', self.listcry)
        print('listcryZ \t\t--> Atomic number of Crystalline species in the substrate:', self.listcryZ)
        print('calibration_params \t--> Calibration parameters:')
        for key, value in self.calibration_params.items():
            print('     {} \t--> {}'.format(key, value))


    # -------------- AUXILIARY ROUTINES 
    # Indices is not used inside the class, but it's useful to write out start.dat (can be fixed better...)
    def indices(self):
        """
        Returns a numpy array containing the indices of elements in a probability matrix for a system with:  
        ncry: number of crystalline species 
        ncov: number of coverage species
        """
        ncry, ncov = len(self.listcry), len(self.listcov)
        x=[0,1,2,3]
        # Find all the possible permutations of elements in x
        import itertools as it
        perm = np.asarray(list(it.product(x, repeat=ncry+ncov)))
        # Exclude those where the sum of all elements is larger than 4 (i.e. the total number of nearest neighbours)
        sumperm = np.sum(perm, axis=1)
        todelete = np.where(sumperm>4)[0]
        # Exclude [0,0,0,0]
        todelete = np.append(todelete, np.where(sumperm==0)[0])
        # Exclude those where the first ncry elements are zero (there is always at least 1 crystal neighbour)
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)==0)[0])
        # Exclude those where the sum of the first ncry elements is larger than 3 (a crystal site cannot have more than 3 crystal neighbours in a tetrahedral configuration) 
        todelete = np.append(todelete, np.where(np.sum(perm[:,:ncry], axis=1)>3)[0])
        # Array of survivors
        inds = np.delete(perm, todelete, axis=0)
        return inds

    def get_dims(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        return np.ones(ncry+ncov, dtype=int)*4    

    # -------------- SETUP EVAPORATION FREQUENCIES
    def get_energetics_evap(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        Etot=np.zeros((ncry,*dims))  # energetics of evaporation in eV
        
        if self.substrate == 'Si':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int)*4
            E0 = np.zeros((ncry,*dims0))
            E0[0,1] = 1*self.calibration_params['Ed2'] # 1coor, no coverage
            E0[0,2] = 2*self.calibration_params['Ed2'] # 2coor, no coverage 
            E0[0,3] = 3*self.calibration_params['Ed2'] # 3coor, no coverage
        elif self.substrate == 'Ge':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int) * 4
            E0 = np.zeros((ncry, *dims0))
            E0[0,1] = 1 * self.calibration_params['Ed2']  # 1coor, no coverage
            E0[0,2] = 2 * self.calibration_params['Ed2']  # 2coor, no coverage
            E0[0,3] = 3 * self.calibration_params['Ed2']  # 3coor, no coverage
        elif self.substrate == 'SiGe02fake':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int) * 4
            E0 = np.zeros((ncry, *dims0))
            E0[0, 1] = 1 * self.calibration_params['Ed2']  # 1coor, no coverage
            E0[0, 2] = 2 * self.calibration_params['Ed2']  # 2coor, no coverage
            E0[0, 3] = 3 * self.calibration_params['Ed2']  # 3coor, no coverage
        elif self.substrate == 'SiGe04fake':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int) * 4
            E0 = np.zeros((ncry, *dims0))
            E0[0,1] = 1 * self.calibration_params['Ed2']  # 1coor, no coverage
            E0[0,2] = 2 * self.calibration_params['Ed2']  # 2coor, no coverage
            E0[0,3] = 3 * self.calibration_params['Ed2']  # 3coor, no coverage
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()

        for i in range(ncry):
            Etot[i] = E0[i]

        return Etot

    def get_PtransE(self):
        
        # Energetics for deposition
        ene = self.get_energetics_evap()
        print('\nEvaporation energetics:\n {}'.format(ene))

        # Prefactor
        P0 = self.calibration_params['P0']
        
        # Calculate Frequencies
        PtransE = P0 * np.exp(-ene/(self.kb_ev*self.temp))
        print('\nEvaporation probability:\n {}'.format(PtransE))

        return PtransE

    # -------------- SETUP DEPOSITION FREQUENCIES
    def get_energetics_depo(self):
        ncry = len(self.listcry)
        Etot=np.zeros((ncry,3))  # energetics of deposition in eV
        
        E0=np.zeros((ncry,3)) 
        if self.substrate == 'Si':
            E0[0,0]= 2*self.calibration_params['Ed1'] # coor 1 
            E0[0,1]= 2*self.calibration_params['Ed2'] # coor 2
            E0[0,2]= 2*self.calibration_params['Ed3'] # coor 3 
        elif self.substrate == 'Ge':
            E0[0,0] = 2 * self.calibration_params['Ed1']  # coor 1
            E0[0,1] = 2 * self.calibration_params['Ed2']  # coor 2
            E0[0,2] = 2 * self.calibration_params['Ed3']  # coor 3
        elif self.substrate == 'SiGe02fake':
            E0[0,0]= 2*self.calibration_params['Ed1'] # coor 1
            E0[0,1]= 2*self.calibration_params['Ed2'] # coor 2
            E0[0,2]= 2*self.calibration_params['Ed3'] # coor 3
        elif self.substrate == 'SiGe04fake':
            E0[0,0]= 2*self.calibration_params['Ed1'] # coor 1
            E0[0,1]= 2*self.calibration_params['Ed2'] # coor 2
            E0[0,2]= 2*self.calibration_params['Ed3'] # coor 3
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
            
        # Set E0 in Etot 
        for i in range(ncry):
            Etot[i] = E0[i]
        
        return Etot

    def dampingf(self, temp):
        A = self.calibration_params['A']
        Tflex = self.calibration_params['Tflex']
        fact = (1 + math.erf((temp - Tflex) / A)) * 0.5  #error function (riscalabile con Tflex e A)
        print('\nDampingfactor:\n {}'.format(fact))
        return fact

    def get_PtransD(self):
        
        # Energetics for deposition
        ene = self.get_energetics_depo()
        print('\nDeposition energetics:\n {}'.format(ene))

        # Prefactor
        P0 = self.calibration_params['P0']
        
        # Calculate Frequencies
        PtransD = P0 * np.exp(-ene/(self.kb_ev*self.Tm)) * self.dampingf(self.temp) #0.5 is not needed here!!!!!
        print('\nDeposition probability:\n {}'.format(PtransD))

        return PtransD




