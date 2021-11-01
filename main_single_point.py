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
import CG_LNA.spectre as sp

#===========================================================================================================================
#------------------------------------ Other Functions ----------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the MOSFET parameters to the circuit_initialization_parameters dictionary
def get_mos_parameters(circuit_initialization_parameters,process_name):
	
	circuit_initialization_parameters['MOS']={}
	circuit_initialization_parameters['MOS']['Process']=process_name
	circuit_initialization_parameters['MOS']['filename']={}
	
	f=open('/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/MOS_Files/'+process_name+'.txt')
	lines=f.readlines()
	f.close()

	# Extracting values from the MOS File
	for i in range(len(lines)):
		line=lines[i][:-1]
		if line=='Vdd':
			circuit_initialization_parameters['MOS']['Vdd']=float(lines[i+1][:-1])
		elif line=='Lmin':
			circuit_initialization_parameters['MOS']['Lmin']=float(lines[i+1][:-1])
		elif line=='u0':
			circuit_initialization_parameters['MOS']['un']=float(lines[i+1][:-1])*1e-4
		elif line=='tox':
			circuit_initialization_parameters['MOS']['tox']=float(lines[i+1][:-1])
		elif line=='vth0':
			circuit_initialization_parameters['MOS']['vt']=float(lines[i+1][:-1])
		elif line=='tt_file':
			circuit_initialization_parameters['MOS']['filename']['tt']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['tt']+=lines[j]
				j+=1
		elif line=='ff_file':
			circuit_initialization_parameters['MOS']['filename']['ff']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['ff']+=lines[j]
				j+=1
		elif line=='ss_file':
			circuit_initialization_parameters['MOS']['filename']['ss']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['ss']+=lines[j]
				j+=1
                
	# Calculating Cox
	eo=8.85*1e-12
	er=3.9
	circuit_initialization_parameters['MOS']['cox']=eo*er/circuit_initialization_parameters['MOS']['tox']

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}
	circuit_initialization_parameters['simulation']['directory']='/home/ee18b028/cadence_project/lna1/'
	circuit_initialization_parameters['simulation']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	circuit_initialization_parameters['simulation']['iip3_type']='basic'		# 'basic' or 'advanced' 

	circuit_initialization_parameters['simulation']['std_temp']=27
	circuit_initialization_parameters['simulation']['pin_fixed']=-65
	circuit_initialization_parameters['simulation']['pin_start']=-70
	circuit_initialization_parameters['simulation']['pin_stop']=-40
	circuit_initialization_parameters['simulation']['pin_points']=6
	circuit_initialization_parameters['simulation']['iip3_calc_points']=3
	circuit_initialization_parameters['simulation']['process_corner']='tt'
	circuit_initialization_parameters['simulation']['conservative']='NO'
	circuit_initialization_parameters['simulation']['w_finger_max']=2e-6

	circuit_initialization_parameters['simulation']['parameters_list']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}

#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}

# ---------- MOSFET Parameters ----------
#get_mos_parameters(circuit_initialization_parameters,'TSMC180')
#get_mos_parameters(circuit_initialization_parameters,'TSMC65')
get_mos_parameters(circuit_initialization_parameters,'TSMC65_2')
#get_mos_parameters(circuit_initialization_parameters,'IBM130')

# ---------- Simulation Conditions ----------
fo=1e9
get_simulation_conditions(circuit_initialization_parameters,fo)

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
	'Rb':274,
	'Rd':260,
	'Io':1270e-6,
	'C1':31.8e-12,
	'C2':124e-12,
	'W':190e-6,
	'Rbias':1000
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE RUN --------------------------------------------
cir=sp.Circuit(circuit_initialization_parameters)
cir.update_circuit(circuit_parameters)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------- OUTPUT PRINT --------------------------------------------
print('Extracted_Parameters\n')
for param_name in cir.extracted_parameters:
	print(param_name,' : ',cir.extracted_parameters[param_name])


#===========================================================================================================================