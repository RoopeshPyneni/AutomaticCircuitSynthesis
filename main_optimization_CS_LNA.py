#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will initialize the optimization_input_parameters and run the complete_optimization file for CS LNA
"""

#===========================================================================================================================
import numpy as np
import complete_optimization as co
import copy


"""
===========================================================================================================================
------------------------------------ OTHER FUNCTIONS ----------------------------------------------------------------------
"""

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
# Function that sets the output conditions to the optimization_input_parameters dictionary
def get_output_conditions(optimization_input_parameters,fo):
	
	optimization_input_parameters['output_conditions']={
		's11_db':-10.0,
		'iip3_dbm':-10.0,
		'gain_db':15.0,
		'gain_delta':1.5,
		'nf_db':1.5,
		'wo':2.0*np.pi*fo,
		'delta_v':0.1,
		'Rs':50,
		'Cload':400e-15,
		's11_db_middle':-15.0
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the circuit_initialization_parameters dictionary
def get_simulation_conditions(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}
	circuit_initialization_parameters['simulation']['standard_parameters']={}

	# Filenames
	circuit_initialization_parameters['simulation']['standard_parameters']['directory']='/home/ee18b028/cadence_project/CS_LNA/'
	circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']='/home/ee18b028/cadence_project/'
	circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']='mos_inductor' # 'ideal', 'series','mos_resistor','mos_capacitor','mos_inductor'
	
	# IIP3 Points
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']='basic'		# 'basic' or 'advanced' 
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=-65
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']=-70
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']=-40
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']=6
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']=3
	circuit_initialization_parameters['simulation']['standard_parameters']['n_harm']=5
	circuit_initialization_parameters['simulation']['standard_parameters']['f_iip3']=1e6

	# Operating frequency points
	circuit_initialization_parameters['simulation']['standard_parameters']['f_list']=[fo-50e6,fo,fo+50e6]
	
	# Other Values
	circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']=27
	circuit_initialization_parameters['simulation']['standard_parameters']['temp_list']=[-40,27,120]
	circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']=['ss','tt','ff']
	circuit_initialization_parameters['simulation']['standard_parameters']['conservative']='NO'
	circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max']=2e-6

	circuit_initialization_parameters['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5,
		'process_corner':'tt'
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the pre_optimization parameters to the optimization_input_parameters dictionary
def get_pre_optimization_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['pre_optimization']={}

	optimization_input_parameters['pre_optimization']['type']='manual'

	optimization_input_parameters['pre_optimization']['I_Rdivider_max']=100e-6
	

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Manual Hand Calculations
	optimization_input_parameters['pre_optimization']['manual_circuit_parameters']={
	'Ld': 9e-9,
	'Cd': 5.489971363501514e-13,
	'W': 0.0025202377594010653,
	'Cg': 1.5458676014581085e-10,
	'Io': 0.0016369962227517215,
	'Rsum': 10000,
	'Rk': 0.95,
	'Ls': 3.182280567225685e-09,
	'Lg': 1.792018017477487e-08,
	'Rb': 5000,
	'Cs': 1.2271773971045355e-11,
	'Wk': 100,
	'Io_k':2.2262529019823066e-06,
	'Wpr':6.951933459876268e-07,
	'wid_Cpr':4.074636060571575e-06,
	}
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Pre Optimization Simulation Parameters
	optimization_input_parameters['pre_optimization']['simulation']={}
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']={}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the optimization parameters to the optimization_input_parameters dictionary
def get_optimization_parameters(optimization_input_parameters,fo,optimization_name):

	# Getting the number of optimization runs
	optimization_input_parameters['optimization']={}
	optimization_input_parameters['optimization']['run']='NO'
	optimization_input_parameters['optimization']['n_runs']=1

	# Getting the type of optimization
	if optimization_name=='LOSS':
		optimization_input_parameters['optimization']['optimization_name']='loss1'
	else:
		optimization_input_parameters['optimization']['optimization_name']='fom1'

	# This dictionary will store the simulation parameters
	optimization_input_parameters['optimization']['simulation']={}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Optimization Run 1
	optimization_input_parameters['optimization'][1]={}
	optimization_input_parameters['optimization'][1]['max_iteration']=25
	optimization_input_parameters['optimization'][1]['alpha_min']=-1
	optimization_input_parameters['optimization'][1]['consec_iter']=-1
	optimization_input_parameters['optimization'][1]['delta_threshold']=0.001
	optimization_input_parameters['optimization'][1]['alpha_mult']=1
	optimization_input_parameters['optimization'][1]['loss_type']=2
	optimization_input_parameters['optimization'][1]['optimization_type']=0
	#optimization_input_parameters['optimization'][1]['optimizing_parameters']=['Lg','Io','W','Ls','Cd','Rk','Cs','Io_k','wid_Cpr','Wpr']
	optimization_input_parameters['optimization'][1]['optimizing_parameters']=['Cd']
	#optimization_input_parameters['optimization'][1]['optimizing_parameters']=['Cd','Wpr','wid_Cpr']
	optimization_input_parameters['optimization'][1]['output_parameters_list']=['Io','gain_db','iip3_dbm','s11_db','Zin_R','Zin_I','nf_db','p_source','gm1','vg1','vd1']
	optimization_input_parameters['optimization'][1]['process_temp_type']='all'

	# NOTES :
	# loss_type will tell how to update the circuit parameters
	#	0 - Consider loss slope of Io if other loss is 0 ; Consider loss slope of others only if other loss of others is not zero
	#	1 - Normal Update
	#	2 - Consider the loss slope of those components whose loss is non-zero

	# optimization_type will tell whether to increase or decrease loss
	# 	0 - we want to reduce the loss
	# 	1 - we want to increase the loss

	optimization_input_parameters['optimization'][1]['loss_weights']={
		'gain_db':1.0/10.0,
		'iip3_dbm':1.0/10.0,
		's11_db':3.0/15.0,
		'nf_db':1.0/2.0,
		'Io':0.0,

		'gain_delta':0.0/10.0,
		'gain_flatness':0.0,
		's11_db_middle':3.0/15.0,
		'gain_delta2':1.0/10.0
	}

	optimization_input_parameters['optimization'][1]['alpha']={}
	optimization_input_parameters['optimization'][1]['alpha']['value']=0.05
	optimization_input_parameters['optimization'][1]['alpha']['type']='Normal' # 'Linear','Log'
	optimization_input_parameters['optimization'][1]['alpha']['start']=0.8
	optimization_input_parameters['optimization'][1]['alpha']['end']=0.05


	# Optimization Simulation Parameters
	optimization_input_parameters['optimization']['simulation'][1]={}
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']={}
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Optimization Run 2
	optimization_input_parameters['optimization'][2]=copy.deepcopy(optimization_input_parameters['optimization'][1])
	optimization_input_parameters['optimization']['simulation'][2]=copy.deepcopy(optimization_input_parameters['optimization']['simulation'][1])
	optimization_input_parameters['optimization'][2]['max_iteration']=5
	optimization_input_parameters['optimization'][2]['optimizing_parameters']=['Lg','Io','W','Ls','Cd','Rk','Cs']
	optimization_input_parameters['optimization'][2]['process_temp_type']='normal'
	
	optimization_input_parameters['optimization'][2]['alpha']['value']=0.005
	
	optimization_input_parameters['optimization'][2]['loss_weights']={
		'gain_db':1.0/10.0,
		'iip3_dbm':1.0/10.0,
		's11_db':3.0/15.0,
		'nf_db':1.0/2.0,
		'Io':100.0,

		'gain_delta':1.0/10.0,
		'gain_flatness':0.0,
		's11_db_middle':3.0/15.0,
		'gain_delta2':0.1/10.0
	}
	
	
#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary
def get_sensitivity_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['sensitivity_analysis']={}
	optimization_input_parameters['sensitivity_analysis']['run']='NO'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Sensitivity Analysis Simulation Parameters
	optimization_input_parameters['sensitivity_analysis']['simulation']={}
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['temp_list']=[-40,27,120]

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary
def get_temperature_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['temperature_analysis']={}

	optimization_input_parameters['temperature_analysis']['run']='NO'

	optimization_input_parameters['temperature_analysis']['start_current']=0.1
	optimization_input_parameters['temperature_analysis']['stop_current']=10
	optimization_input_parameters['temperature_analysis']['n_current']=1

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['temperature_analysis']['simulation']={}
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['temp_list']=[-40,0,40,80,120]

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the process analysis parameters to the optimization_input_parameters dictionary
def get_process_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['process_analysis']={}
	optimization_input_parameters['process_analysis']['run']='NO'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['process_analysis']['simulation']={}
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['temp_list']=[-40,0,40,80,120]
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['process_corner']=['ss','tt','ff']

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the frequency analysis parameters to the optimization_input_parameters dictionary
def get_frequency_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['frequency_analysis']={}

	optimization_input_parameters['frequency_analysis']['run']='YES'

	optimization_input_parameters['frequency_analysis']['plot_type']='linear' # 'log'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Frequency Analysis Simulation Parameters
	optimization_input_parameters['frequency_analysis']['simulation']={}
	optimization_input_parameters['frequency_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['frequency_analysis']['simulation']['standard_parameters']['f_list']=[0.8e9,0.9e9,1.0e9,1.1e9,1.2e9]


#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary
def get_iip3_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['iip3_analysis']={}
	optimization_input_parameters['iip3_analysis']['run']='NO'

	optimization_input_parameters['iip3_analysis']['pin_start']=-70
	optimization_input_parameters['iip3_analysis']['pin_stop']=-40
	optimization_input_parameters['iip3_analysis']['n_pin']=16
	optimization_input_parameters['iip3_analysis']['n_points']=5
	optimization_input_parameters['iip3_analysis']['freq_array']=[fo]

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['iip3_analysis']['simulation']={}
	optimization_input_parameters['iip3_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['iip3_analysis']['simulation']['standard_parameters']['f_list']=[0.8e9,0.9e9,1.0e9,1.1e9,1.2e9]

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the frequency analysis parameters to the optimization_input_parameters dictionary
def get_circuit_parameter_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['circuit_parameter_analysis']={}

	optimization_input_parameters['circuit_parameter_analysis']['run']='YES'
	optimization_input_parameters['circuit_parameter_analysis']['n_runs']=1
	
	optimization_input_parameters['circuit_parameter_analysis'][0]={}
	optimization_input_parameters['circuit_parameter_analysis'][0]['parameter_name']='Cd'
	optimization_input_parameters['circuit_parameter_analysis'][0]['parameter_select_type']='relative'
	optimization_input_parameters['circuit_parameter_analysis'][0]['start']=0.8
	optimization_input_parameters['circuit_parameter_analysis'][0]['stop']=1.2
	optimization_input_parameters['circuit_parameter_analysis'][0]['n_value']=11
	optimization_input_parameters['circuit_parameter_analysis'][0]['sweep_type']='linear' # 'log'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Frequency Analysis Simulation Parameters
	optimization_input_parameters['circuit_parameter_analysis']['simulation']={}
	optimization_input_parameters['circuit_parameter_analysis']['simulation']['standard_parameters']={}
	#optimization_input_parameters['circuit_parameter_analysis']['simulation']['standard_parameters']['f_list']=[0.8e9,0.9e9,1.0e9,1.1e9,1.2e9]


"""
===========================================================================================================================
------------------------------------ MAIN PROGRAM -------------------------------------------------------------------------
"""

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}
optimization_input_parameters={}
optimization_name='LOSS'

"""
~~~~~~~~~~ MAIN PARAMETERS ~~~~~~~~~
"""

# ---------- MOSFET Parameters ----------
#get_mos_parameters(circuit_initialization_parameters,'TSMC180')
get_mos_parameters(circuit_initialization_parameters,'TSMC65')
#get_mos_parameters(circuit_initialization_parameters,'IBM130')

# ---------- Output Conditions ----------
fo=1e9
get_output_conditions(optimization_input_parameters,fo)

# ---------- Simulation Conditions ----------
get_simulation_conditions(circuit_initialization_parameters,fo)

# ---------- Pre Optimization Parameters ----------
get_pre_optimization_parameters(optimization_input_parameters,fo)

# ---------- Optimization Parameters ----------
get_optimization_parameters(optimization_input_parameters,fo,optimization_name)

"""
~~~~~~~~~~ ANALYSIS PARAMETERS ~~~~~~~~~
"""

# ---------- Sensitivity Analysis Parameters ----------
get_sensitivity_analysis_parameters(optimization_input_parameters,fo)

# ---------- Temperature Analysis Parameters ----------
get_temperature_analysis_parameters(optimization_input_parameters,fo)

# ---------- Process Analysis Parameters ----------
get_process_analysis_parameters(optimization_input_parameters,fo)

# ---------- IIP3 Analysis Parameters ----------
get_iip3_analysis_parameters(optimization_input_parameters,fo)

# ---------- Frequency Analysis Parameters ----------
get_frequency_analysis_parameters(optimization_input_parameters,fo)

# ---------- Circuit Parameter Analysis Parameters ----------
get_circuit_parameter_analysis_parameters(optimization_input_parameters,fo)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE NAMES ------------------------------------------

optimization_input_parameters['filename']={}
optimization_input_parameters['filename']['run_status']='/home/ee18b028/Optimization/Simulation_Results/run_status.txt'

f_directory='/home/ee18b028/Optimization/Simulation_Results/CS_LNA/'


file_choose='S' # 'S' to run a single time; 'M' to run multiple times

optimization_input_parameters['optimization']['run']='YES' #'YES'
optimization_input_parameters['temperature_analysis']['run']='NO'
optimization_input_parameters['sensitivity_analysis']['run']='NO'
optimization_input_parameters['process_analysis']['run']='NO'
optimization_input_parameters['iip3_analysis']['run']='NO'
optimization_input_parameters['frequency_analysis']['run']='NO' #'YES'
optimization_input_parameters['circuit_parameter_analysis']['run']='NO' #'YES'

if file_choose=='S':

	# ------- Set Any Additional Parameters Here --------
	filename=f_directory+'4_Optimization_1'						# SET THE FILENAME HERE
	# ------- Set Any Additional Parameters Here --------
	

	# ------- DON'T CHANGE THESE LINES -------------
	optimization_input_parameters['filename']['output']=filename
	co.complete_optimization(circuit_initialization_parameters,optimization_input_parameters,'CS_LNA')			
	# ------- DON'T CHANGE THESE LINES -------------	


if file_choose=='M':
	for i in range(0,5):	# SET NUMBER OF ITERATIONS HERE

		# ------- Set Any Additional Parameters Here --------
		filename=f_directory+'tsmc_str_1'+str(i)			# SET THE FILENAME HERE
		fo=1e9
		wo=2*np.pi*fo
		optimization_input_parameters['output_conditions']['wo']=wo
		iip_mtd=['basic_hb','basic_pss','advanced_hb','advanced_pss','hb_sweep']
		optimization_input_parameters['simulation']['iip3_method']=iip_mtd[i]
		# ------- Set Any Additional Parameters Here --------



		# ------- DON'T CHANGE THESE LINES -------------
		optimization_input_parameters['filename']['output']=filename
		co.complete_optimization(optimization_input_parameters)
		# ------- DON'T CHANGE THESE LINES -------------

#===========================================================================================================================
