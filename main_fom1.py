#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Main Code:
"""
#===========================================================================================================================
import numpy as np
import complete_optimization as co
import data_plot as dp


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
#-------------------------------------- OUTPUT CONDITIONS --------------------------------------

fo=1e9
optimization_input_parameters['output_conditions']={
	's11_db':-15.0,
	'iip3_dbm':-5.0,
	'gain_db':10.0,
	'nf_db':4.0,
	'wo':2.0*np.pi*fo,
	'delta_v':0.1,
	'Rs':50
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------------------------- SPECTRE SIMULATOR --------------------------------------

optimization_input_parameters['simulation']={}
optimization_input_parameters['simulation']['directory']='/home/ee18b028/cadence_project/lna1/'
optimization_input_parameters['simulation']['iip3_method']='hb_sweep' # Options : 'basic_hb', 'basic_pss', 'advanced_hb', 'advanced_pss', 'hb_sweep'
optimization_input_parameters['simulation']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
optimization_input_parameters['simulation']['iip3_calc_points']=5

optimization_input_parameters['simulation']['parameters_list']={
	'pin':-65,
	'fund_2':fo+1e6,
	'fund_1':fo,
	'cir_temp':27,
	'n_harm':15,
	'pin_start':-80,
	'pin_stop':-60,
	'n_pin':10
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------- Parameters for PRE OPTIMIZATION -------------------------------

optimization_input_parameters['pre_optimization']={}

optimization_input_parameters['pre_optimization']['Step1b_Limit']=5
optimization_input_parameters['pre_optimization']['Step2_Limit']=5
optimization_input_parameters['pre_optimization']['vdsat_reqd']=0.09

optimization_input_parameters['pre_optimization']['type']=1
optimization_input_parameters['pre_optimization']['gmrs_threshold']=0.2
optimization_input_parameters['pre_optimization']['vdsat_threshold']=0.02

optimization_input_parameters['pre_optimization']['C1_threshold']=100
optimization_input_parameters['pre_optimization']['C2_threshold']=100
optimization_input_parameters['pre_optimization']['Rbias_threshold']=100

#~~~~~~~~~~~~~~~~~~~~~~~~~
# Manual Hand Calculations
optimization_input_parameters['pre_optimization']['manual_circuit_parameters']={
	'Rb':614,
	'Rd':316,
	'Io':900e-6,
	'C1':318e-12,
	'C2':106e-12,
	'W':61.5e-6,
	'Rbias':148
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#--------------------------------- Parameters for OPTIMIZATION ---------------------------------

optimization_input_parameters['optimization']={}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters for Optimization
optimization_input_parameters['optimization']['max_iteration']	=300
optimization_input_parameters['optimization']['alpha_min']	=-1
optimization_input_parameters['optimization']['consec_iter']	=-1

optimization_input_parameters['optimization']['delta_threshold']=0.001
optimization_input_parameters['optimization']['alpha_mult']	=1
optimization_input_parameters['optimization']['loss_type']	=1
optimization_input_parameters['optimization']['update_check']	=0

optimization_input_parameters['optimization']['optimizing_parameters']=['Rb','Rd','Io','W']
optimization_input_parameters['optimization']['output_parameters_list']=['Io','gain_db','iip3_dbm','s11_db','nf_db','p_source','gm1','vdsat','vg','vd','vs']

optimization_input_parameters['optimization']['optimization_name']='loss1'
optimization_input_parameters['optimization']['optimization_type']=0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assigning values to the loss weights
loss_weights={}
loss_weights['gain_db']=1/10.0	
loss_weights['iip3_dbm']=1/5.0	
loss_weights['s11_db']=1/15.0	
loss_weights['nf_db']=1/4.0	
loss_weights['Io']=1000	
optimization_input_parameters['optimization']['loss_weights']=loss_weights


#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assigning Values of Alpha
alpha_parameters={}
alpha_parameters['common']=0.2
alpha_parameters['Rb']=1
alpha_parameters['Rd']=1
alpha_parameters['W']=1
alpha_parameters['Io']=1
optimization_input_parameters['optimization']['alpha']={}
optimization_input_parameters['optimization']['alpha']['values']=alpha_parameters

optimization_input_parameters['optimization']['alpha']['type']='Normal'
optimization_input_parameters['optimization']['alpha']['start']=0.8
optimization_input_parameters['optimization']['alpha']['end']=0.05



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditions for acceptable solution
optimization_input_parameters['acceptable_solution']={}
optimization_input_parameters['acceptable_solution']['s11_db']=-15
optimization_input_parameters['acceptable_solution']['gain_db']=6
optimization_input_parameters['acceptable_solution']['iip3_dbm']=-15
optimization_input_parameters['acceptable_solution']['nf_db']=10
optimization_input_parameters['acceptable_solution']['p_source']=10e-3



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----------------------------------------- FILE NAMES ------------------------------------------

optimization_input_parameters['filename']={}
optimization_input_parameters['filename']['run_status']='/home/ee18b028/Optimization/Simulation_Results/run_status.txt'

f_directory='/home/ee18b028/Eldo_files/Opt_Code/Simulation_Results/FOM/'

file_choose='M' # 'S' to run a single time; 'M' to run multiple times


if file_choose=='S':
	# ------- Set Any Additional Parameters Here --------
	filename=f_directory+'hc_2_optimization'						# SET THE FILENAME HERE
	optimization_input_parameters['max_iteration']=500
	# ------- Set Any Additional Parameters Here --------
	

	# ------- DON'T CHANGE THESE LINES -------------
	optimization_input_parameters['filename']['output']=filename
	co.complete_optimization(optimization_input_parameters)
	dp.plot_complete(optimization_input_parameters)
	# ------- DON'T CHANGE THESE LINES -------------

if file_choose=='M':
	for i in range(1,2):
		# ------- Set Any Additional Parameters Here --------
		filename=f_directory+'fom_newHC_ibm_run_'+str(i)		# SET THE FILENAME HERE
		fo=1e9*i
		wo=2*np.pi*fo
		optimization_input_parameters['output_conditions']['wo']=wo
		# ------- Set Any Additional Parameters Here --------

		
		# ------- DON'T CHANGE THESE LINES -------------
		optimization_input_parameters['filename']['output']=filename
		co.complete_optimization(optimization_input_parameters)
		dp.plot_complete(optimization_input_parameters)
		# ------- DON'T CHANGE THESE LINES -------------


#===========================================================================================================================