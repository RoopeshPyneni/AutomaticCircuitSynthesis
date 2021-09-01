#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Main Code:
"""
#===========================================================================================================================
import numpy as np
import math
import spectre as sp


#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating a dictionary with the optimization parameters
optimization_input_parameters={}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Parameters for MOSFET SELECTION -------------------------------

optimization_input_parameters['MOS']={}

#"""
optimization_input_parameters['MOS']['filename']='/home/ee18b028/cadence_project/tsmc018.scs'
optimization_input_parameters['MOS']['Type']='NMOS'
optimization_input_parameters['MOS']['Lmin']=0.18*1e-6
optimization_input_parameters['MOS']['Vdd']=1.8
#"""

"""
optimization_input_parameters['MOS']['filename']='/home/ee18b028/cadence_project/ibm013.scs'
optimization_input_parameters['MOS']['Type']='NMOS'
optimization_input_parameters['MOS']['Lmin']=0.13*1e-6
optimization_input_parameters['MOS']['Vdd']=1.3
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------------------------- SPECTRE SIMULATOR --------------------------------------

optimization_input_parameters['simulation']={}
optimization_input_parameters['simulation']['directory']='/home/ee18b028/cadence_project/lna1/'
optimization_input_parameters['simulation']['basic_circuit']='basic_parameters'
optimization_input_parameters['simulation']['iip3_circuit']='iip3_hb'
optimization_input_parameters['simulation']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
optimization_input_parameters['simulation']['iip3_type']='advanced'		# 'basic' or 'advanced' 

optimization_input_parameters['simulation']['std_temp']=27
optimization_input_parameters['simulation']['pin_fixed']=-65
optimization_input_parameters['simulation']['pin_start']=-70
optimization_input_parameters['simulation']['pin_stop']=-40
optimization_input_parameters['simulation']['pin_points']=16
optimization_input_parameters['simulation']['iip3_calc_points']=5

fo=1e9
optimization_input_parameters['simulation']['parameters_list']={
	'pin':-65,
	'fund_2':fo+1e6,
	'fund_1':fo,
	'cir_temp':27,
	'n_harm':15
}

optimization_input_parameters['simulation']['cir_writing_dict']={
	'wid':'W',
	'cur0':'Io',
	'Resb':'Rb',
	'Resd':'Rd',
	'cap1':'C1',
	'cap2':'C2',
	'Resbias':'Rbias'
}


circuit_parameters={
	'Rb':262,
	'Rd':1310,
	'Io':437e-6,
	'C1':159e-12,
	'C2':1960e-12,
	'W':184e-6,
	'Rbias':8.1
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE RUN --------------------------------------------

sp.write_MOS_parameters(optimization_input_parameters)
extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------- OUTPUT PRINT --------------------------------------------
print('Extracted_Parameters\n')
for param_name in extracted_parameters:
	print(param_name,' : ',extracted_parameters[param_name])



#===========================================================================================================================

