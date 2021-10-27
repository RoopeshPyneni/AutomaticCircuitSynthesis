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
	optimization_input_parameters['MOS']['Process']=process_name
	optimization_input_parameters['MOS']['filename']={}
	
	f=open('/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/MOS_Files/'+process_name+'.txt')
	lines=f.readlines()
	f.close()

	# Extracting values from the MOS File
	for i in range(len(lines)):
		line=lines[i][:-1]
		if line=='Vdd':
			optimization_input_parameters['MOS']['Vdd']=float(lines[i+1][:-1])
		elif line=='Lmin':
			optimization_input_parameters['MOS']['Lmin']=float(lines[i+1][:-1])
		elif line=='u0':
			optimization_input_parameters['MOS']['un']=float(lines[i+1][:-1])*1e-4
		elif line=='tox':
			optimization_input_parameters['MOS']['tox']=float(lines[i+1][:-1])
		elif line=='vth0':
			optimization_input_parameters['MOS']['vt']=float(lines[i+1][:-1])
		elif line=='tt_file':
			optimization_input_parameters['MOS']['filename']['tt']=''
			j=i+1
			while lines[j][:-1]!='':
				optimization_input_parameters['MOS']['filename']['tt']+=lines[j]
				j+=1
		elif line=='ff_file':
			optimization_input_parameters['MOS']['filename']['ff']=''
			j=i+1
			while lines[j][:-1]!='':
				optimization_input_parameters['MOS']['filename']['ff']+=lines[j]
				j+=1
		elif line=='ss_file':
			optimization_input_parameters['MOS']['filename']['ss']=''
			j=i+1
			while lines[j][:-1]!='':
				optimization_input_parameters['MOS']['filename']['ss']+=lines[j]
				j+=1
                
	# Calculating Cox
	eo=8.85*1e-12
	er=3.9
	optimization_input_parameters['MOS']['cox']=eo*er/optimization_input_parameters['MOS']['tox']

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions(optimization_input_parameters,fo):
	
	optimization_input_parameters['simulation']={}
	optimization_input_parameters['simulation']['directory']='/home/ee18b028/cadence_project/lna1/'
	optimization_input_parameters['simulation']['basic_circuit']='s11_parameters_tsmc_65'
	optimization_input_parameters['simulation']['iip3_circuit']='iip3_hb_tsmc_65_test'
	optimization_input_parameters['simulation']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	optimization_input_parameters['simulation']['iip3_type']='basic'		# 'basic' or 'advanced' 

	optimization_input_parameters['simulation']['std_temp']=27
	optimization_input_parameters['simulation']['pin_fixed']=-65
	optimization_input_parameters['simulation']['pin_start']=-70
	optimization_input_parameters['simulation']['pin_stop']=-40
	optimization_input_parameters['simulation']['pin_points']=16
	optimization_input_parameters['simulation']['iip3_calc_points']=5
	optimization_input_parameters['simulation']['process_corner']='tt'
	optimization_input_parameters['simulation']['conservative']='NO'
	optimization_input_parameters['simulation']['w_finger_max']=2e-6

	optimization_input_parameters['simulation']['parameters_list']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':10
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
#get_mos_parameters(optimization_input_parameters,'TSMC180')
#get_mos_parameters(optimization_input_parameters,'TSMC65')
get_mos_parameters(optimization_input_parameters,'TSMC65_2')
#get_mos_parameters(optimization_input_parameters,'IBM130')

# ---------- Simulation Conditions ----------
fo=1e9
get_simulation_conditions(optimization_input_parameters,fo)

"""
circuit_parameters={
	'Rb':272,
	'Rd':360,
	'Io':666e-6,
	'C1':31.8e-12,
	'C2':163e-12,
	'W':275e-6,
	'Rbias':1000
}
"""

circuit_parameters={
	'Rb':273,
	'Rd':265,
	'Io':1300e-6,
	'C1':31.8e-12,
	'C2':126e-12,
	'W':184e-6,
	'Rbias':1000
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
	
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- CHANGING TO NEW SCHEMATIC -------------------------------------
	
optimization_input_parameters['simulation']['basic_circuit']='basic_parameters_tsmc_65'
optimization_input_parameters['simulation']['iip3_circuit']='iip3_hb_tsmc_65'
#circuit_parameters['C1']=31.8e-12

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE RUN --------------------------------------------

sp.write_MOS_parameters(optimization_input_parameters)
extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------- OUTPUT PRINT --------------------------------------------
print('Extracted_Parameters\n')
for param_name in extracted_parameters:
	print(param_name,' : ',extracted_parameters[param_name])

"""



#===========================================================================================================================
