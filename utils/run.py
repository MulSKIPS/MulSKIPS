import os
import numpy as np
import pymulskips as pm #  pymulskips should be in your $PYTHONPATH

"""
This is a simple script to run PVD growth of a flat SiC-3C surface using the pymulskips library
"""

execpath = '../mulskips-source/'
runpath = os.getenv('PWD')+'/SiC3C-surface/'

# Data extracted from the excel file from Peter Wellmann "Temperature vs Growth rate.xlsx" 25 November 2020
# New sight glass. # Power in kW, T in C
data_samples = {}
data_samples['151'] = {'Power': 5.45, 'T_seed_middle': 1986.26, 'T_seed_edge': 0, 'T_source_middle': 1994.14, 'T_source_edge': 1994.14, 'Exp-Growth-Rate': 101}
data_samples['148'] = {'Power': 5.85, 'T_seed_middle': 2040.19, 'T_seed_edge': 0, 'T_source_middle': 2048.12, 'T_source_edge': 2048.12, 'Exp-Growth-Rate': 208}
data_samples['149'] = {'Power': 5.90, 'T_seed_middle': 2046.70, 'T_seed_edge': 0, 'T_source_middle': 2054.63, 'T_source_edge': 2054.63, 'Exp-Growth-Rate': 230}
data_samples['152'] = {'Power': 6.00, 'T_seed_middle': 2059.56, 'T_seed_edge': 0, 'T_source_middle': 2067.51, 'T_source_edge': 2067.51, 'Exp-Growth-Rate': 231}
data_samples['153'] = {'Power': 6.10, 'T_seed_middle': 2072.24, 'T_seed_edge': 0, 'T_source_middle': 2080.20, 'T_source_edge': 2080.20, 'Exp-Growth-Rate': 287}

# Select process ID
key = '151'

# Get temperatures
T_zero_kelvin = 273.15
Tsource = (data_samples[key]['T_source_middle'] + data_samples[key]['T_source_edge']) * 0.5 + T_zero_kelvin # [K]
Tseed_center = data_samples[key]['T_seed_middle'] + T_zero_kelvin # [K]

# Setup PVD process
p = pm.process.PVD(substrate='SiC-3C', precursors=['Si', 'Si2C', 'SiC2'], 
	calibration_type='avrov', Tsource=Tsource, Tseed_center=Tseed_center)

# Run MulSKIPS
pm.setuprun.run_mulskips(execpath, runpath, Simulation='F', mp=p, 
	Seed_box=[120, 360, 360], PtransZig=1.0, OutMolMol=100000, IterMax=1000000)
