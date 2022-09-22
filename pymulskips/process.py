import time
import math
import numpy as np
import os,shutil,subprocess,sys, copy
import cantera as ct


# ----------------------------------------------------------------------------------------------------------------------
# CLASS FOR PVD PROCESSES
# ----------------------------------------------------------------------------------------------------------------------

class PVD():

    def __init__(self, 
        substrate='SiC-3C', precursors=['Si', 'Si2C', 'SiC2'], 
        calibration_type='avrov', Tsource=1500, Tseed_center=1560,
        calibration_params=None, miller=100):
        
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

        self.miller = miller

        # By default here there is no coverage
        self.listcov = []
        self.listcovZ = []

        # Some constants related to the gas phase
        self.Rgas = ct.gas_constant * 1e-3       # Gas constant J/mol/K
        self.NA = ct.avogadro * 1e-3       # Avogadro Number [mol^-1]
        self.kb_ev = ct.boltzmann / ct.electron_charge     # [eV/K]
        self.h = ct.planck / ct.electron_charge # eV*s

        # Molar masses of elements (Si, C, H, Cl, B, P, Ge...)
        self.M_Si = ct.Element('Si').weight*1e-3  #Kg/mol
        self.M_H = ct.Element('H').weight*1e-3  #Kg/mol
        self.M_C = ct.Element('C').weight*1e-3 #Kg/mol
        self.M_Cl = ct.Element('Cl').weight*1e-3 #Kg/mol
        self.M_Ge = ct.Element('Ge').weight*1e-3 #Kg/mol
        self.M_B = ct.Element('B').weight*1e-3 #Kg/mol
        self.M_P = ct.Element('P').weight*1e-3 #Kg/mol


        # Some constants related to the solid growing substrate
        if self.substrate == 'SiC-3C': 
            self.listcry = ['Si', 'C']
            self.mass = self.M_Si + self.M_C # [kg/mol]
            self.rho = 3210 # kg/m^3    
            self.a0 = pow(4*self.mass/self.NA/self.rho,1.0/3.0) # [m] side of the cubic conventional cell of 3C-SiC containing four Si-C dimers
            self.KMC_sf = 4.63/12.0   # KMC Super-Lattice parameter (angstrom)
        elif self.substrate == 'Si': 
            self.listcry = ['Si']
            self.mass = self.M_Si # [kg/mol]
            self.rho = 2329  # [kg/m^3]
            self.a0 = pow(8 * self.mass / self.NA / self.rho,1.0 / 3.0)  # [m]  #lattice constant of Si (diamond structure, containing eight atoms per cell)
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)
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
        # Atomic numbers
        self.listcryZ = [ct.Element(el).atomic_number for el in self.listcry]
        # Below is an approiximated rho_surf. For Miller indices differientiation, we should include them for each face in the lines above 
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

    def __init__(self, expdata, calibration_params=None, miller=100, gas=None):
        
        print('Reading input from provided dictionary')
        self.expdata = expdata
        self.substrate = expdata['substrate'] 
        self.precursors = expdata['precursors'] 
        self.temp = expdata['temperature']
        self.partial_pressures = expdata['partial_pressures']

        self.miller = miller

        self.gas = gas

        # Some constants related to the gas phase
        self.Rgas = ct.gas_constant * 1e-3       # Gas constant J/mol/K
        self.NA = ct.avogadro * 1e-3       # Avogadro Number [mol^-1]
        self.kb_ev = ct.boltzmann / ct.electron_charge     # [eV/K]
        self.h = ct.planck / ct.electron_charge # eV*s

        # Molar masses of elements (Si, C, H, Cl, B, P, Ge...)
        self.M_Si = ct.Element('Si').weight*1e-3  #Kg/mol
        self.M_H = ct.Element('H').weight*1e-3  #Kg/mol
        self.M_C = ct.Element('C').weight*1e-3 #Kg/mol
        self.M_Cl = ct.Element('Cl').weight*1e-3 #Kg/mol
        self.M_Ge = ct.Element('Ge').weight*1e-3 #Kg/mol
        self.M_B = ct.Element('B').weight*1e-3 #Kg/mol
        self.M_P = ct.Element('P').weight*1e-3 #Kg/mol

        # Some constants related to the solid growing substrate
        if self.substrate == 'SiC-3C': 
            self.listcry = ['Si', 'C']
            self.mass = self.M_Si + self.M_C # [kg/mol]
            self.rho = 3210 # kg/m^3    
            self.a0 = pow(4*self.mass/self.NA/self.rho,1.0/3.0) # [m] side of the cubic conventional cell of 3C-SiC containing four Si-C dimers
            self.KMC_sf = 4.63/12.0   # KMC Super-Lattice parameter (angstrom)
        elif self.substrate == 'Si': 
            self.listcry = ['Si']
            self.mass = self.M_Si # [kg/mol]
            self.rho = 2329  # [kg/m^3]
            self.a0 = pow(8 * self.mass / self.NA / self.rho,1.0 / 3.0)  # [m]  #lattice constant of Si (diamond structure, containing eight atoms per cell)
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)
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
        # Atomic numbers
        self.listcryZ = [ct.Element(el).atomic_number for el in self.listcry]
        # Below is an approiximated rho_surf. For Miller indices differientiation, we should include them for each face in the lines above 
        self.rho_surf = self.rho * self.a0 # Kg/m^2


        # Define process substrate + precursors    
        if self.gas:
            print('Using cantera equilibrium major species as precursors:\n', self.precursors)
            self.precursor_masses = {}
            tmp = []
            for mol in self.precursors: 
                tmp.extend(list(gas.species(mol).composition.keys()))
                self.precursor_masses[mol] = gas.molecular_weights[gas.species_index(mol)]*1e-3 # kg/mol
            self.listcov = list(set(tmp)-set(self.listcry))
            self.listcovZ = [ct.Element(el).atomic_number for el in self.listcov]

            calibration_params_keys = ['Edep', 'Eev', 'Eabs', 'Edes', 
                'deltaE', 'kevap_A', 'kevap_E', 'scalef', 'scalefcov']


        else:
            if self.substrate == 'Si':
                if self.precursors == ['SiH4', 'H2']:
                    # Some constants related to the gas-phase molecules
                    M_SiH4 = self.M_Si + 4*self.M_H # [kg/mol]
                    M_H2 = 2*self.M_H # [kg/mol]
                    self.precursor_masses = {'SiH4': M_SiH4, 'H2': M_H2}
                    self.listcov = ['H']
                    self.listcovZ = [1]

                    calibration_params_keys = ['Edep', 'Eev', 'Eabs', 'Edes', 'deltaH', 'deltaCl', 
                        'Ereact_SiH4', 'Ereact_H2', 'k_Si_SiH4', 'k_H_SiH4', 'k_H_H2', 'scalef', 
                        'aads', 'bads', 's0', 'Ades', 'Aev', 'Bev', 'Adep']

                else:
                    print('ERROR: CVD process for substrate {} with precursors {} is not implemented'.format(self.substrate, self.precursors))
                    sys.exit()
            else:
                print('ERROR: CVD process for substrate {} is not implemented'.format(self.substrate))
                sys.exit()

        # Check calibration params
        if calibration_params is None:
            print('ERROR: Please provide calibration parameters dictionary...')
            sys.exit()

        for key in calibration_params_keys:
            if key not in list(calibration_params.keys()):
                print('  ERROR: please provide calibration parameter {} '.format(key))
                sys.exit()

        # for key, value in calibration_params.items():
        #     if key not in calibration_params_keys:
        #         print('  ERROR: please provide calibration parameter {} '.format(key))
        #         sys.exit()
        # Assign calibration params
        self.calibration_params = copy.deepcopy(calibration_params)


        print('*** Initializing CVD process for {} substrate with precursors {}'.format(self.substrate, self.precursors))
        print('mass \t\t\t--> Substrate Mass [Kg/mol]:', self.mass)
        print('rho \t\t\t--> Substrate Density [Kg/m^3]:', self.rho)
        print('rho_surf \t\t--> Substrate Surface density [Kg/m^2]:', self.rho_surf)
        print('at_density \t\t--> Substrate Surface atomic density [at/m^2]:', self.at_density)
        print('temp \t\t\t--> Substrate temperature [K]:', self.temp)
        print('precursor_masses \t--> Precursor Masses [kg/mol]:\n', self.precursor_masses)
        print('partial_pressures \t\t--> Precursors partial pressures [Pa]:\n', self.partial_pressures)
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
        #     p_dep = self.partial_pressures[mol]   # Pascal
        #     coeff = pow(2.*math.pi*self.precursor_masses[mol]*self.Rgas*self.temp,-0.5)
        p_dep = self.partial_pressures[mol]   # Pascal
        coeff = pow(2.*math.pi*self.precursor_masses[mol]*self.Rgas*self.temp,-0.5)
        flux = p_dep * coeff # mol/s/m^2 
        return flux

    def get_fluxes(self):
        fluxes = {}
        for mol in self.precursors:
            fluxes[mol] = self.get_flux(mol) # mol/s/m^2
        return fluxes



    # -------------- SETUP DEPOSITION FREQUENCIES

    def s0(self, mol): 
        return self.calibration_params['alpha'][mol]/(1+self.calibration_params['kd0kr0ratio'][mol] \
                * np.exp(-(self.calibration_params['EdminusEr'][mol])/(self.kb_ev*self.temp))) #non-dimensional

    def get_fluxes_depo(self):
        
        ncry, ncov = len(self.listcry), len(self.listcov)

        # Get precursors fluxes on surface
        loc_fluxes = self.get_fluxes() # mol/s/m^2
        
        # Get atomic fluxes on surface
        tot_fluxes = {}
        if self.gas: # use cantera results for precursors
            for el in self.listcry+self.listcov:
                tot_fluxes[el] = 0
                for mol in self.precursors:
                    if el in self.gas.species(mol).composition:
                        stoy = self.gas.species(mol).composition[el]
                        tot_fluxes[el] += loc_fluxes[mol] * stoy * self.s0(mol)  # /self.NA   # Na lo abbiamo messo per compensare il fattore in alpha (18Feb2022) 
                        # tot_fluxes[el] += loc_fluxes[mol] * stoy * self.s0(mol)  # * ( 1 - np.sum(theta) / self.at_density)     # POSSIBILE STRATEGIA ALTERNATIVA

        else:
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
                    tot_fluxes['Si'] = loc_fluxes['SiH4'] * self.calibration_params['k_Si_SiH4'] * reactf_SiH4  # Silicon
                    tot_fluxes['H'] = loc_fluxes['SiH4'] * self.calibration_params['k_H_SiH4'] * reactf_SiH4 + \
                                    loc_fluxes['H2'] * self.calibration_params['k_H_H2'] * reactf_H2  # Hydrogen
                else:
                    print('ERROR: not implemented'); sys.exit()
            else:
                print('ERROR: Substrate not implemented'); sys.exit()

        return tot_fluxes

    def get_energetics_depo(self):
        ncry = len(self.listcry)
        Etot=np.zeros((ncry,3))  # energetics of deposition in eV
        for iel, el in enumerate(self.listcry):
            Etot[iel] = self.calibration_params['Edep'][el] # list for coor 1, 2 and 3           
        return Etot

    def get_PtransD(self):
        
        # Energetics for deposition
        ene = self.get_energetics_depo()
        print('\nDeposition energetics:\n {}'.format(ene))

        # Prefactor
        #freq_0 = self.mass / self.rho_surf # m^2 / mol
        # freq_0 = self.NA / self.at_density  # m^2/mol

        # Fluxes
        fluxes_depo = self.get_fluxes_depo()

        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        PtransDtot = np.zeros((ncry,3))
        Bolt = np.zeros((ncry,3))
        for i in range(ncry):
            # ll = np.zeros(ncry,dtype=int)
            # if ncry==1:
                # ll[-(1+i)]=1
            # elif ncry>1:
                # ll[-(1+i)]=2
            # enecoo2 = ene[tuple(np.r_[i,ll])]  # select element in matrix corresponding to coor=2; for ncry=2, consider that the 2 neighbours should be of the "other" species
            Bolt[i]=np.exp(-ene[i]/(self.kb_ev * self.temp)) # Make sure Bolt(n=2) = 1 --> Edep2 = 0
            PtransDtot[i] = self.calibration_params['scalef'] * fluxes_depo[self.listcry[i]] * Bolt[i]  # at/s
            print(fluxes_depo[self.listcry[i]], self.listcry[i])   

        return PtransDtot


    # -------------- SETUP EVAPORATION FREQUENCIES

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
        If ncov is not zero, return deltaEtot: energy barrier for evaporation in eV, induced by coverage
        Note:
        H --> increases deltaEtot
        Cl --> decreases deltaEtot
        """
        # Optionally provide an extra adjustment 
        ncry, ncov = len(self.listcry), len(self.listcov)
        
        # Check that all keys in deltaE belong to listcov
        if not all(item in self.listcov for item in list(self.calibration_params['deltaE'].keys())):
            print('ERROR: Coverage species {} is not implemented'.format(self.listcov[ii]))         
            sys.exit()

        # # Check that H increases deltaEtot (deltaE values are all >=0)
        # if not all(np.array(self.calibration_params['deltaE']['H']) >= 0):
        #     print('ERROR: Please provide deltaE values >= 0 for H')         
        #     sys.exit()
        # # Check that Cl decreases deltaEtot (deltaE values are all <=0)
        # if not all(np.array(self.calibration_params['deltaE']['Cl']) <= 0):
        #     print('ERROR: Please provide deltaE values <= 0 for Cl')         
        #     sys.exit()

        deltaEtot=0
        if ncov!=0:
            covn = indici[ncry:]
            for ii,nn in enumerate(covn):
                # deltaEtot += self.calibration_params['deltaE'][self.listcov[ii]] * nn  # # deltaE should be a dict with scalar items
                if nn!=0: deltaEtot += self.calibration_params['deltaE'][self.listcov[ii]][nn-1] # deltaE here should be a dict with dim(3)array items  
        return deltaEtot

    def get_dims(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        return np.ones(ncry+ncov, dtype=int)*4    

    def get_energetics_evap(self):
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        Etot=np.zeros((ncry,*dims))  # energetics of evaporation in eV
        
        for iel, el in enumerate(self.listcry):
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int)*4
            E0 = np.zeros((ncry,*dims0))
            E0[0,1] = self.calibration_params['Eev'][el][0] # 1coor, no coverage
            E0[0,2] = self.calibration_params['Eev'][el][1] # 2coor, no coverage 
            E0[0,3] = self.calibration_params['Eev'][el][2] # 3coor, no coverage
            
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

    def get_kevap(self):
        # Get atomic fluxes on surface
        tot_kevap = {}
        for el in self.listcry+self.listcov:
            tot_kevap[el] = 0
            for mol in self.precursors:
                if el in self.gas.species(mol).composition:
                    stoy = self.gas.species(mol).composition[el]
                    tot_kevap[el] += stoy * self.calibration_params['kevap_A'][mol] * np.exp(-self.calibration_params['kevap_E'][mol]/(self.kb_ev*self.temp))
        return tot_kevap

    def get_PtransE(self):    

        # Energetics for evaporation
        ene = self.get_energetics_evap()
        print('\nEvaporation energetics:\n {}'.format(ene))

        # Evaporation rates
        kevap = self.get_kevap()

        # Calculate Frequencies
        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        PtransEtot = np.zeros((ncry,*dims))
        Bolt = np.zeros((ncry,*dims))
        for i in range(ncry):
            Bolt[i] = np.exp(-ene[i]/ (self.kb_ev * self.temp))
            PtransEtot[i] = self.calibration_params['scalef'] * kevap[self.listcry[i]] * Bolt[i] # at/s
            print(kevap[self.listcry[i]], self.listcry[i])

        return PtransEtot



    # -------------- SETUP ABSORPTION FREQUENCIES

    def get_energetics_abs(self):
        ncov = len(self.listcov)
        Etot=np.zeros((ncov,3))  # energetics of absorption in eV
        for iel, el in enumerate(self.listcov):
            Etot[iel] = self.calibration_params['Eabs'][el] # list for coor 1, 2 and 3   # <---- TO BE CALIBRATED        
        return Etot

    def get_PtransAbs(self):    
        # Initialize probability matrix with dimension depending on NCrystal and NCov
        # Energetics for absorption
        ene = self.get_energetics_abs()
        print('\nAbsorption energetics:\n {}'.format(ene))

        # Fluxes
        fluxes_depo = self.get_fluxes_depo()

        # Calculate Frequencies
        ncov = len(self.listcov)
        PtransAbstot = np.zeros((ncov,3))
        Bolt = np.zeros((ncov,3))
        for i in range(ncov):
            Bolt = np.exp(-ene/ (self.kb_ev * self.temp))
            PtransAbstot[i] = self.calibration_params['scalefcov'] * fluxes_depo[self.listcov[i]] * Bolt[i] # at/s
            print(fluxes_depo[self.listcov[i]], self.listcov[i])
        return PtransAbstot 

    # -------------- SETUP DESORPTION FREQUENCIES

    def get_energetics_des(self):
        ncov = len(self.listcov)
        Etot=np.zeros((ncov,3))  # energetics of absorption in eV
        for iel, el in enumerate(self.listcov):
            Etot[iel] = self.calibration_params['Edes'][el] # list for coor 1, 2 and 3   # <---- TO BE CALIBRATED
        return Etot

    def get_PtransDes(self):    
        # Initialize probability matrix with dimension depending on NCrystal and NCov

        # Energetics for absorption
        ene = self.get_energetics_des()
        print('\nDesorption energetics:\n {}'.format(ene))

        # Evaporation rates
        kevap = self.get_kevap()

        ncov = len(self.listcov)
        PtransDestot = np.zeros((ncov,3))
        Bolt = np.zeros((ncov,3))
        for i in range(ncov):
            Bolt[i]=np.exp(-ene[i]/(self.kb_ev * self.temp))   # non normalizziamo ocn enecoo2 qui, perché quel ragionamento sul coor=2 che riproduce l'andamento medio vale quando il prefattore è funzione del flusso
            PtransDestot[i] = self.calibration_params['scalefcov'] * kevap[self.listcov[i]] * Bolt[i] # at/s
            print(kevap[self.listcov[i]], self.listcov[i])
        return PtransDestot





# ----------------------------------------------------------------------------------------------------------------------
# CLASS FOR LA PROCESSES
# ----------------------------------------------------------------------------------------------------------------------

class PhaseChange(): # when we have Ge and SiGE calibration, rename this as melting_MulskipsProcess()

    def __init__(self, material, temp=None, calibration_params=None):  
        
        self.substrate = material
        if temp is not None:
            self.temp = temp
        if calibration_params is not None:
            self.calibration_params = calibration_params

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
            self.Tm = [1688] # K
            self.KMC_sf = 5.43/12.0   # KMC Super-Lattice parameter (angstrom)

            if calibration_params is None: #values used in npj article #mytag
                calibration_params = {
                    # Energy variation with coordination
                    'delta': 0.03, # [eV]
                    # Energy barrier for deposition of solid Si in a site with coor 2
                    'Ed2': 0.96, # [eV]
                    # Prefactor
                    'P0': 1.33e17,
                    # constant in damping exp for depo probs
                    'A': 280,
                    'Tflex': 1080, # [K]
                    # Molar fraction of species
                    'X0':1.0,
                    'X':1.0
                }

        elif self.substrate == 'Ge':
            self.mass = 0.07263 # [kg/mol]
            self.rho = 5323  # [kg/m^3]
            self.listcry = ['Ge']
            self.listcryZ = [32]
            self.Tm = [1210]  # [K] #modified as list #mytag
            self.KMC_sf = 5.66/12.0  # KMC Super-Lattice parameter (angstrom)

            if calibration_params is None: #calibration to be confirmed!!!!!!!!!!!!!!!!!!!!!!!!! #mytag
                calibration_params = {
                    # Energy variation with coordination
                    'delta': 0.03*0.65, # [eV]
                    # Energy barrier for deposition of solid Si in a site with coor 2
                    'Ed2': 0.65, # [eV]
                    # Prefactor
                    'P0': 2.596e16,
                    # constant in damping exp for depo probs
                    'A': 210,
                    'Tflex': 900, # [K]
                    # Molar fraction of species
                    'X0':1.0,
                    'X':1.0
                }
        elif self.substrate == 'SiGe_X':  # interpolations strictly valid for X in [0.0,0.4] and X=1
            self.mass = 0.07263 * self.calibration_params['X0'][1] + 0.0280855 * (1.0 - self.calibration_params['X0'][1])  # [kg/mol]
            self.rho = 2329 + 3493 * self.calibration_params['X0'][1] - 499 * pow(self.calibration_params['X0'][1],2)  # [kg/m^3]
            self.listcry = ['Si','Ge']
            self.listcryZ = [14, 32]
            self.Tm = [1688, 1210]
            #self.Tm = 1210 * self.calibration_params['X'][1] + 1688 * (1.0 - self.calibration_params['X'][1])  # K
            #self.KMC_sf = (5.43 * (1.0 - self.calibration_params['X'][1]) + 5.66 * self.calibration_params['X'][1]) / 12.0  # KMC Super-Lattice parameter (angstrom)
            self.KMC_sf = 5.43/12.0  # KMC Super-Lattice parameter (angstrom)
            # NB: KMC_sf for strained SiGe is similar to Si...we need it to be fixed to Si anyways for LA where X is variable! 
            # Then for melt depth we can use a weighted average depending on x (in turn evaluated as a func. of xyz at the end of the process) 

        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()

        a0 = pow(8*self.mass/self.NA/self.rho,1.0/3.0) # [m] lattice parameter for crystalline Si (diamond structure, 8 atoms per cell)
        print('lattice parameter for substrate:', a0)
        self.rho_surf = self.rho * a0 # Kg/m^2
                
        # Check
        calibration_params_keys = ['Ed2', 'delta', 'P0', 'A', 'Tflex', 'X', 'X0'] #mytag_2903
        for key, value in calibration_params.items():
            if key not in calibration_params_keys:
                # print('  ERROR: please provide calibration parameter {} '.format(key))
                # sys.exit()
                print('  WARNING: unknown calibration parameter {} will not be taken into account'.format(key))
        for key in calibration_params_keys:
            if key not in calibration_params:
                print('  ERROR: please provide calibration parameter {} '.format(key))
                sys.exit()

        # calibration_params_keys = ['Ed1', 'Ed2', 'Ed3', 'P0', 'A', 'Tflex']
        # for key, value in calibration_params.items():
        #     if key not in calibration_params_keys:
        #         print('  ERROR: please provide calibration parameter {} '.format(key))
        #         sys.exit()
        # self.calibration_params = calibration_params

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
            E0[0,1] = 1*self.calibration_params['Ed2'][0] # 1coor, no coverage
            E0[0,2] = 2*self.calibration_params['Ed2'][0] # 2coor, no coverage
            E0[0,3] = 3*self.calibration_params['Ed2'][0] # 3coor, no coverage
        elif self.substrate == 'Ge':
            # Energetics for evaporation without coverage neighbours
            dims0 = np.ones(ncry, dtype=int) * 4
            E0 = np.zeros((ncry, *dims0))
            E0[0,1] = 1 * self.calibration_params['Ed2'][0]  # 1coor, no coverage
            E0[0,2] = 2 * self.calibration_params['Ed2'][0] # 2coor, no coverage
            E0[0,3] = 3 * self.calibration_params['Ed2'][0]  # 3coor, no coverage
        elif self.substrate == 'SiGe_X':  # ncry=2 behaving like Si(1-x)+Ge(x)
            # Energetics for evaporation without coverage neighbours
            # da modificare in funzione di self.calibration_params['X']
            dims0 = np.ones(ncry, dtype=int) * 4
            E0 = np.zeros((ncry, *dims0))
            alpha = 0.063 #0.06 #0.07 calib2 #0.06 new calib
            beta = 0.085 #0.08  #0.2 calib2 #0.08 new calib
            # ---------------------------------------------------------------------------------------------------------------
            Ed2 = self.calibration_params['Ed2']
            E0[0, 1, 0] = 1 * Ed2[0]                                    # Si, 1Si NN, no coverage 
            E0[0, 0, 1] = (1-alpha) * 1 * (Ed2[0] + Ed2[1]) / 2         # Si, 1Ge NN, no coverage
            E0[0, 1, 1] = (1-alpha) * 2 * (3 * Ed2[0] + Ed2[1]) / 4     # Si, 1Si 1Ge NN, no coverage
            E0[0, 2, 0] = 2 * Ed2[0]                                    # Si, 2Si NN, no coverage 
            E0[0, 0, 2] = (1-alpha) * 2 * (Ed2[0] + Ed2[1]) / 2         # Si, 2Ge NN, no coverage 
            E0[0, 2, 1] = (1-alpha) * 3 * (5 * Ed2[0] + Ed2[1]) / 6     # Si, 2Si 1Ge NN, no coverage
            E0[0, 1, 2] = (1-alpha) * 3 * (2 * Ed2[0] + Ed2[1]) / 3     # Si, 1Si 2Ge NN, no coverage
            E0[0, 3, 0] = 3 * Ed2[0]                                    # Si, 3Si NN, no coverage 
            E0[0, 0, 3] = (1-alpha) * 3 * (Ed2[0] + Ed2[1]) / 2         # Si, 3Ge NN, no coverage
            E0[1, 1, 0] = (1+2*beta) * 1 * (Ed2[1] + Ed2[0]) / 2          # Ge, 1Si NN, no coverage
            E0[1, 0, 1] = 1 * Ed2[1]                                    # Ge, 1Ge NN, no coverage 
            E0[1, 1, 1] = (1+beta) * 2 * (3 * Ed2[1] + Ed2[0]) / 4      # Ge, 1Si 1Ge NN, no coverage 
            E0[1, 2, 0] = (1+2.5*beta) * 2 * (Ed2[1] + Ed2[0]) / 2          # Ge, 2Si NN, no coverage
            E0[1, 0, 2] = 2 * Ed2[1]                                    # Ge, 2Ge NN, no coverage 
            E0[1, 2, 1] = (1+beta) * 3 * (2 * Ed2[1] + Ed2[0]) / 3      # Ge, 2Si 1Ge NN, no coverage
            E0[1, 1, 2] = (1+beta) * 3 * (5 * Ed2[1] + Ed2[0]) / 6      # Ge, 1Si 2Ge NN, no coverage
            E0[1, 3, 0] = (1+4*beta) * 3 * (Ed2[1] + Ed2[0]) / 2          # Ge, 3Si NN, no coverage
            E0[1, 0, 3] = 3 * Ed2[1]                                    # Ge, 3Ge NN, no coverage 
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()

        for i in range(ncry):
            Etot[i] = E0[i]

        print('\nEvaporation energetics:\n {}'.format(Etot))

        return Etot

    def get_PtransE(self):
        
        # Energetics for deposition
        ene = self.get_energetics_evap()
        #print('\nEvaporation energetics:\n {}'.format(ene))

        # Prefactor #mytag
        #P0 = self.calibration_params['P0']
        #PtransE = P0 * np.exp(-ene/(self.kb_ev*self.temp))

        ncry, ncov = len(self.listcry), len(self.listcov)
        dims = self.get_dims()
        P0 = np.zeros((ncry, *dims))
        #P0 = np.zeros((len(self.listcry), 4, 1)) #mytag (si dovrebbe poter cancellare!!!!!!!)

        for i, j in enumerate(self.listcry):
            P0[i] = self.calibration_params['P0'][i]

        # Calculate Frequencies
        PtransE = np.multiply(P0, np.exp(-ene / (self.kb_ev * self.temp)))
        print('\nEvaporation probability:\n {}'.format(PtransE))

        return PtransE

    # -------------- SETUP DEPOSITION FREQUENCIES
    def get_energetics_depo(self):
        ncry = len(self.listcry)
        Etot=np.zeros((ncry,3))  # energetics of deposition in eV
        
        E0=np.zeros((ncry,3)) 
        if self.substrate == 'Si':
            E0[0,0] = 2 * (self.calibration_params['Ed2'][0] - self.calibration_params['delta'][0])  # coor 1
            E0[0,1] = 2 * self.calibration_params['Ed2'][0]  # coor 2
            E0[0,2] = 2 * (self.calibration_params['Ed2'][0] + self.calibration_params['delta'][0])  # coor 3
        elif self.substrate == 'Ge':
            E0[0,0] = 2 * (self.calibration_params['Ed2'][0] - self.calibration_params['delta'][0])  # coor 1
            E0[0,1] = 2 * self.calibration_params['Ed2'][0]  # coor 2
            E0[0,2] = 2 * (self.calibration_params['Ed2'][0] + self.calibration_params['delta'][0])  # coor 3
        elif self.substrate == 'SiGe_X': # ncry=2
            Ed2 = self.calibration_params['Ed2']
            delta = self.calibration_params['delta']
            E0[0,0] = 2 * (Ed2[0] - delta[0])     # Si coor 1
            E0[0,1] = 2 * Ed2[0]                  # Si coor 2
            E0[0,2] = 2 * (Ed2[0] + delta[0])     # Si coor 3
            E0[1,0] = 2 * (Ed2[1] - delta[1])     # Ge coor 1
            E0[1,1] = 2 * Ed2[1]                  # Ge coor 2
            E0[1,2] = 2 * (Ed2[1] + delta[1])     # Ge coor 3 
        else:
            print('ERROR: Substrate {} is not implemented'.format(self.substrate))
            sys.exit()
            
        # Set E0 in Etot
        for i in range(ncry):
            Etot[i] = E0[i]

        print('\nDeposition energetics:\n {}'.format(Etot))
        
        return Etot


    def dampingf(self, temp): #mytag
        A = np.zeros(len(self.listcry))
        Tflex = np.zeros(len(self.listcry))
        fact = np.zeros((len(self.listcry), 3))
        for i, Z in enumerate(self.listcry):
            A[i] = self.calibration_params['A'][i]
            Tflex[i] = self.calibration_params['Tflex'][i]
            fact[i] = (1 + math.erf((temp - Tflex[i]) / A[i])) * 0.5  # error function (riscalabile con Tflex e A)
        return fact

    def get_PtransD(self): #mytag
        
        # Energetics for deposition
        ene = self.get_energetics_depo()
        #print('\nDeposition energetics:\n {}'.format(ene))

        # Prefactor
        P0 = np.zeros((len(self.listcry), 1))
        Tm = np.zeros((len(self.listcry), 1))
        for i, Z in enumerate(self.listcry):
            P0[i] = self.calibration_params['P0'][i] * self.calibration_params['X'][i]
            Tm[i] = self.Tm[i]
            # Calculate Frequencies
            PtransD = P0 * np.exp(-ene / (self.kb_ev * Tm)) * self.dampingf(self.temp)
            print('\ni:',i)
            print('Tm(i):', Tm[i])
            print('P0(i):', P0[i])
            print('Dampingfactor(i):', self.dampingf(self.temp)[i])

        return PtransD
