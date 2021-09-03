#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: main_single_point.py
File Description 	: This file will run spectre and extract the parameters for a single point

Functions structure in this file:
	--> get_mos_parameters
	--> get_simulation_conditions
	
"""

#===========================================================================================================================
import spectre as sp

#===========================================================================================================================
#------------------------------------ Other Functions ----------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the MOSFET parameters to the optimization_input_parameters dictionary
def get_mos_parameters(optimization_input_parameters,process_name):
	
	optimization_input_parameters['MOS']={}
	
	if process_name=='TSMC018':
		optimization_input_parameters['MOS']['filename']='/home/ee18b028/cadence_project/tsmc018.scs'
		optimization_input_parameters['MOS']['Type']='NMOS'
		optimization_input_parameters['MOS']['Lmin']=0.18*1e-6
		optimization_input_parameters['MOS']['Vdd']=1.8
	else:
		optimization_input_parameters['MOS']['filename']='/home/ee18b028/cadence_project/ibm013.scs'
		optimization_input_parameters['MOS']['Type']='NMOS'
		optimization_input_parameters['MOS']['Lmin']=0.13*1e-6
		optimization_input_parameters['MOS']['Vdd']=1.3

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions(optimization_input_parameters,fo):
	
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
	optimization_input_parameters['simulation']['pin_points']=6
	optimization_input_parameters['simulation']['iip3_calc_points']=3

	optimization_input_parameters['simulation']['parameters_list']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
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

#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

# Creating a dictionary with the optimization parameters
optimization_input_parameters={}

# ---------- MOSFET Parameters ----------
get_mos_parameters(optimization_input_parameters,'TSMC018')
#get_mos_parameters(optimization_input_parameters,'IBM013')

# ---------- Simulation Conditions ----------
fo=1e9
get_simulation_conditions(optimization_input_parameters,fo)


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