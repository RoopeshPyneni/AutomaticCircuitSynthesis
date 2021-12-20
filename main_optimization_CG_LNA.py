#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will initialize the optimization_input_parameters and run the complete_optimization file
"""

#===========================================================================================================================
import numpy as np
import complete_optimization as co

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
# Function that sets the output conditions to the optimization_input_parameters dictionary
def get_output_conditions(optimization_input_parameters,fo):
	
	optimization_input_parameters['output_conditions']={
		's11_db':-15.0,
		'iip3_dbm':-5.0,
		'gain_db':10.0,
		'nf_db':4.0,
		'wo':2.0*np.pi*fo,
		'delta_v':0.1,
		'Rs':50
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the circuit_initialization_parameters dictionary
def get_simulation_conditions(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}
	circuit_initialization_parameters['simulation']['standard_parameters']={}

	circuit_initialization_parameters['simulation']['standard_parameters']['directory']='/home/ee18b028/cadence_project/lna1/'
	circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']='basic'		# 'basic' or 'advanced' 

	circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']=27
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=-65
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']=-70
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']=-40
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']=6
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']=3
	circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']='tt'
	circuit_initialization_parameters['simulation']['standard_parameters']['conservative']='NO'
	circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max']=2e-6

	circuit_initialization_parameters['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the pre_optimization parameters to the optimization_input_parameters dictionary
def get_pre_optimization_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['pre_optimization']={}

	optimization_input_parameters['pre_optimization']['Step1b_Limit']=5
	optimization_input_parameters['pre_optimization']['Step2_Limit']=5
	optimization_input_parameters['pre_optimization']['vdsat_reqd']=0.07

	optimization_input_parameters['pre_optimization']['type']='manual'
	optimization_input_parameters['pre_optimization']['gmrs_threshold']=0.2
	optimization_input_parameters['pre_optimization']['vdsat_threshold']=0.02

	optimization_input_parameters['pre_optimization']['C1_threshold']=100
	optimization_input_parameters['pre_optimization']['C2_threshold']=100
	optimization_input_parameters['pre_optimization']['Rbias_threshold']=100
	optimization_input_parameters['pre_optimization']['Rbias_minimum']=1000

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Manual Hand Calculations
	optimization_input_parameters['pre_optimization']['manual_circuit_parameters']={
		'Rb':274,
		'Rd':360,
		'Io':666e-6,
		'C1':31.8e-12,
		'C2':163e-12,
		'W':275e-6,
		'Rbias':1000
	}
	

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Pre Optimization Simulation Parameters
	optimization_input_parameters['pre_optimization']['simulation']={}
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']={}

	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['std_temp']=27
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_start']=-70
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['pin_points']=6
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['iip3_calc_points']=3
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['pre_optimization']['simulation']['standard_parameters']['conservative']='NO'

	optimization_input_parameters['pre_optimization']['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the optimization parameters to the optimization_input_parameters dictionary
def get_optimization_parameters(optimization_input_parameters,fo,optimization_name):

	optimization_input_parameters['optimization']={}
	optimization_input_parameters['optimization']['run']='YES'
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
	optimization_input_parameters['optimization'][1]['max_iteration']=600
	optimization_input_parameters['optimization'][1]['alpha_min']=-1
	optimization_input_parameters['optimization'][1]['consec_iter']=-1
	optimization_input_parameters['optimization'][1]['delta_threshold']=0.001
	optimization_input_parameters['optimization'][1]['alpha_mult']=1
	optimization_input_parameters['optimization'][1]['loss_type']=0
	optimization_input_parameters['optimization'][1]['update_check']=0
	optimization_input_parameters['optimization'][1]['optimization_type']=0
	optimization_input_parameters['optimization'][1]['optimizing_parameters']=['Rb','Rd','Io','W','C1','C2']
	optimization_input_parameters['optimization'][1]['output_parameters_list']=['Io','gain_db','iip3_dbm','s11_db','s12_db','s21_db','s22_db','k','nf_db','p_source','gm1','vdsat','vg','vd','vs']

	optimization_input_parameters['optimization'][1]['loss_weights']={
		'gain_db':1/10.0,	
		'iip3_dbm':1/5.0,	
		's11_db':1/15.0,
		'nf_db':1/4.0,
		'Io':1000
	}

	optimization_input_parameters['optimization'][1]['alpha']={}
	#optimization_input_parameters['optimization'][1]['alpha']['values']={
	#	'common':0.05,
	#	'Ld':1,
	#	'Lg':1,
	#	'Ls':1,
	#	'W':1,
	#	'Io':1
	#}
	optimization_input_parameters['optimization'][1]['alpha']['value']=0.2
	optimization_input_parameters['optimization'][1]['alpha']['type']='Normal'
	optimization_input_parameters['optimization'][1]['alpha']['start']=0.8
	optimization_input_parameters['optimization'][1]['alpha']['end']=0.05


	# Optimization Simulation Parameters
	optimization_input_parameters['optimization']['simulation'][1]={}
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']={}
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['std_temp']=27
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_start']=-70
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['pin_points']=16
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['iip3_calc_points']=5
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['optimization']['simulation'][1]['standard_parameters']['conservative']='NO'

	optimization_input_parameters['optimization']['simulation'][1]['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}

	"""
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Optimization Run 2
	optimization_input_parameters['optimization'][2]={}
	optimization_input_parameters['optimization'][2]['max_iteration']=600
	optimization_input_parameters['optimization'][2]['alpha_min']=-1
	optimization_input_parameters['optimization'][2]['consec_iter']=-1
	optimization_input_parameters['optimization'][2]['delta_threshold']=0.001
	optimization_input_parameters['optimization'][2]['alpha_mult']=1
	optimization_input_parameters['optimization'][2]['loss_type']=0
	optimization_input_parameters['optimization'][2]['update_check']=0
	optimization_input_parameters['optimization'][2]['optimization_type']=0
	optimization_input_parameters['optimization'][2]['optimizing_parameters']=['Rb','Rd','Io','W','C1','C2']
	optimization_input_parameters['optimization'][2]['output_parameters_list']=['Io','gain_db','iip3_dbm','s11_db','s12_db','s21_db','s22_db','k','nf_db','p_source','gm1','vdsat','vg','vd','vs']

	optimization_input_parameters['optimization'][2]['loss_weights']={
		'gain_db':1/10.0,	
		'iip3_dbm':1/5.0,	
		's11_db':1/15.0,
		'nf_db':1/4.0,
		'Io':1000
	}

	optimization_input_parameters['optimization'][2]['alpha']={}
	#optimization_input_parameters['optimization'][2]['alpha']['values']={
	#	'common':0.05,
	#	'Ld':1,
	#	'Lg':1,
	#	'Ls':1,
	#	'W':1,
	#	'Io':1
	#}
	optimization_input_parameters['optimization'][2]['alpha']['value']=0.2
	optimization_input_parameters['optimization'][2]['alpha']['type']='Normal'
	optimization_input_parameters['optimization'][2]['alpha']['start']=0.8
	optimization_input_parameters['optimization'][2]['alpha']['end']=0.05


	# Optimization Simulation Parameters
	optimization_input_parameters['optimization']['simulation'][2]={}
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']={}
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['std_temp']=27
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_start']=-70
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['pin_points']=16
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['iip3_calc_points']=5
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['optimization']['simulation'][2]['standard_parameters']['conservative']='NO'

	optimization_input_parameters['optimization']['simulation'][2]['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}
	"""
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Conditions for acceptable solution
	if optimization_name=='FOM':
		optimization_input_parameters['acceptable_solution']={}
		optimization_input_parameters['acceptable_solution']['s11_db']=-15
		optimization_input_parameters['acceptable_solution']['gain_db']=6
		optimization_input_parameters['acceptable_solution']['iip3_dbm']=-15
		optimization_input_parameters['acceptable_solution']['nf_db']=10
		optimization_input_parameters['acceptable_solution']['p_source']=10e-3

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary
def get_sensitivity_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['sensitivity_analysis']={}
	optimization_input_parameters['sensitivity_analysis']['run']='YES'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['sensitivity_analysis']['simulation']={}
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['std_temp']=27
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['pin_start']=-70
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['pin_points']=16
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['iip3_calc_points']=5
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['sensitivity_analysis']['simulation']['standard_parameters']['conservative']='YES'

	optimization_input_parameters['sensitivity_analysis']['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':15
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the temperature analysis parameters to the optimization_input_parameters dictionary
def get_temperature_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['temperature_analysis']={}

	optimization_input_parameters['temperature_analysis']['run']='YES'

	optimization_input_parameters['temperature_analysis']['start_temp']=-40
	optimization_input_parameters['temperature_analysis']['stop_temp']=120
	optimization_input_parameters['temperature_analysis']['n_temp']=5

	optimization_input_parameters['temperature_analysis']['start_current']=0.1
	optimization_input_parameters['temperature_analysis']['stop_current']=10
	optimization_input_parameters['temperature_analysis']['n_current']=1


	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['temperature_analysis']['simulation']={}
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['std_temp']=27
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_start']=-70
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['pin_points']=16
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['iip3_calc_points']=5
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['temperature_analysis']['simulation']['standard_parameters']['conservative']='YES'

	optimization_input_parameters['temperature_analysis']['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':15
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the process analysis parameters to the optimization_input_parameters dictionary
def get_process_analysis_parameters(optimization_input_parameters,fo):

	optimization_input_parameters['process_analysis']={}
	optimization_input_parameters['process_analysis']['run']='YES'

	#~~~~~~~~~~~~~~~~~~~~~~~~~
	# Temperature Analysis Simulation Parameters
	optimization_input_parameters['process_analysis']['simulation']={}
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']={}

	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_type']='basic'
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['std_temp']=27
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_fixed']=-65
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_start']=-70
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_stop']=-40
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['pin_points']=16
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['iip3_calc_points']=5
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['process_corner']='tt'
	optimization_input_parameters['process_analysis']['simulation']['standard_parameters']['conservative']='YES'

	optimization_input_parameters['process_analysis']['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':15
	}


#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}
optimization_input_parameters={}
optimization_name='LOSS'

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

# ---------- Sensitivity Analysis Parameters ----------
get_sensitivity_analysis_parameters(optimization_input_parameters,fo)

# ---------- Temperature Analysis Parameters ----------
get_temperature_analysis_parameters(optimization_input_parameters,fo)

# ---------- Process Analysis Parameters ----------
get_process_analysis_parameters(optimization_input_parameters,fo)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE NAMES ------------------------------------------

optimization_input_parameters['filename']={}
optimization_input_parameters['filename']['run_status']='/home/ee18b028/Optimization/Simulation_Results/run_status.txt'

f_directory='/home/ee18b028/Optimization/Simulation_Results/'+str(optimization_name)+'/'


file_choose='S' # 'S' to run a single time; 'M' to run multiple times


if file_choose=='S':

	# ------- Set Any Additional Parameters Here --------
	filename=f_directory+'Test_Simulation_Parameters'						# SET THE FILENAME HERE
	# ------- Set Any Additional Parameters Here --------
	

	# ------- DON'T CHANGE THESE LINES -------------
	optimization_input_parameters['filename']['output']=filename
	co.complete_optimization(circuit_initialization_parameters,optimization_input_parameters,'CG_LNA')			
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
