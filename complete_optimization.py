#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: complete_optimization.py
File Description 	: This file will perform all the optimization steps by calling functions in different files

Functions structure in this file:
	--> save_input_results_initial
	--> save_output_results_initial
	--> save_time_results
	--> complete_optimization

"""

#===========================================================================================================================
import datetime
import optimization as op
import common_functions as cf
import CG_LNA.pre_optimization as pr1
import CS_LNA.pre_optimization as pr2
import temperature_analysis as ta
import process_analysis as pa
import CG_LNA.spectre as sp1
import CS_LNA.spectre as sp2
import os

#===========================================================================================================================
#------------------------------------ File Writing Functions ---------------------------------------------------------------

#-----------------------------------------------------------------
# Function that stores input data of the simulation ( output conditions, MOS Parameters, Filenames, Simulation Conditions )
# Inputs  : optimization_input_parameters
# Outputs : NONE
def save_input_results_initial(optimization_input_parameters):

	# Creating the folder path
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	# Opening the filename
	filename=filename+str('/input_data.txt')
	f=open(filename,'w')

	# Saving Filenames
	f.write('\n\n---------------------- Filenames -----------------------')
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))

	"""
	# Saving MOS Parameters
	f.write('\n\n---------------------- MOS Parameters -----------------------')
	f.write('\nProcess Name	:'+str(optimization_input_parameters['MOS']['Process']))
	f.write('\nVdd      	:'+str(optimization_input_parameters['MOS']['Vdd']))
	f.write('\nLmin     	:'+str(optimization_input_parameters['MOS']['Lmin']))
	f.write('\nun	    	:'+str(optimization_input_parameters['MOS']['un']))
	f.write('\ntox  	   	:'+str(optimization_input_parameters['MOS']['tox']))
	f.write('\ncox	     	:'+str(optimization_input_parameters['MOS']['cox']))
	f.write('\nvth0	     	:'+str(optimization_input_parameters['MOS']['vt']))
	"""

	# Saving Output Conditions
	f.write('\n\n---------------------- Output Conditions -----------------------')
	for name in optimization_input_parameters['output_conditions']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['output_conditions'][name],3))

	"""
	# Saving Simulation Results
	f.write('\n\n---------------------- Simulation Conditions -----------------------')
	f.write('\nDirectory      :'+str(optimization_input_parameters['simulation']['directory']))
	f.write('\nBasic Filename :'+str(optimization_input_parameters['simulation']['basic_circuit']))
	f.write('\nIIP3 Filename  :'+str(optimization_input_parameters['simulation']['iip3_circuit']))
	f.write('\nTCSH Filename  :'+str(optimization_input_parameters['simulation']['tcsh']))
	f.write('\nStandard Temp  :'+str(optimization_input_parameters['simulation']['std_temp']))
	f.write('\nPin Fixed      :'+str(optimization_input_parameters['simulation']['pin_fixed']))
	f.write('\nPin Start      :'+str(optimization_input_parameters['simulation']['pin_start']))
	f.write('\nPin Stop       :'+str(optimization_input_parameters['simulation']['pin_stop']))
	f.write('\nPin Points     :'+str(optimization_input_parameters['simulation']['pin_points']))
	f.write('\nIIP3 Calculation Points :'+str(optimization_input_parameters['simulation']['iip3_calc_points']))
	
	# Parameter List for simulation
	for name in optimization_input_parameters['simulation']['parameters_list']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['simulation']['parameters_list'][name],3))
	"""
	
	f.close()

#-----------------------------------------------------------------
# Function that creates a file that stores output data of the simulation
# Inputs  : optimization_input_parameters
# Outputs : NONE
def save_output_results_initial(optimization_input_parameters):
	
	# Creating the folder path
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	# Opening the filename
	filename=filename+str('/output_data.txt')
	f=open(filename,'w')

	# Storing the Filenames
	f.write('\n\n---------------------- Filenames -----------------------')
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))
	f.write('\n\n\n')

	f.close()

#-----------------------------------------------------------------
# Function that stores the time results
# Inputs  : timing_results, optimization_input_parameters
# Outputs : NONE
def save_time_results(timing_results,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	
	filename=filename+str('/output_data.txt')
	f=open(filename,'a')

	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Timing Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	for optimization_name in timing_results:
		f.write('\n\n---------- '+str(optimization_name)+' ----------')
		
		if 'start' in timing_results[optimization_name]:
			f.write('\nStart    : '+str(timing_results[optimization_name]['start']))
			f.write('\nEnd      : '+str(timing_results[optimization_name]['stop']))
			f.write('\nDuration : '+str(timing_results[optimization_name]['stop']-timing_results[optimization_name]['start']))
		
		else:
			for run_number in timing_results[optimization_name]:
				f.write('\n\nRun Number : '+str(run_number))
				f.write('\nStart    : '+str(timing_results[optimization_name][run_number]['start']))
				f.write('\nEnd      : '+str(timing_results[optimization_name][run_number]['stop']))
				f.write('\nDuration : '+str(timing_results[optimization_name][run_number]['stop']-timing_results[optimization_name][run_number]['start']))
	
	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the MOS File Calculations
# Inputs  : mos_parameters, circuit_initialization_parameters
# Outputs : NONE
def save_mos_results(mos_parameters,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOS Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	for param_name in mos_parameters:
		f.write('\n'+str(param_name)+': '+cf.num_trunc(mos_parameters[param_name],3))
	
	f.close()

#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that performs the complete optimization process
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def complete_optimization(circuit_initialization_parameters,optimization_input_parameters,name):

	# Calculating Starting Time
	timing_results={}
	timing_results['complete_analysis']={}
	timing_results['complete_analysis']['start']=datetime.datetime.now()

	# Saving the optimization input and output results initially
	save_input_results_initial(optimization_input_parameters)
	save_output_results_initial(optimization_input_parameters)

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'w')
	f.write('Filename : '+optimization_input_parameters['filename']['output']+'\n\n')
	f.close()

	#======================================================== MOSFET PARAMETERS ==================================================================================================

	# Writing the MOSFET File Location to .scs file
	if name=='CG_LNA':
		cir=sp1.Circuit(circuit_initialization_parameters)
	else:
		cir=sp2.Circuit(circuit_initialization_parameters)
	save_mos_results(cir.mos_parameters,optimization_input_parameters)

	#======================================================== PRE OPTIMIZATION ===================================================================================================

	if name=='CG_LNA':
		pr1.pre_optimization(cir,optimization_input_parameters,timing_results)
	else:
		pr2.pre_optimization(cir,optimization_input_parameters,timing_results)

	#======================================================== OPTIMIZATION =======================================================================================================

	op.main_opt(cir,optimization_input_parameters,timing_results)

	#======================================================== TEMPERATURE ANALYSIS ===============================================================================================

	ta.temperature_analysis(cir,optimization_input_parameters,timing_results)

	#========================================================== PROCESS ANALYSIS =================================================================================================

	pa.process_analysis(cir,optimization_input_parameters,timing_results)
	
	#======================================================== AFTER OPTIMIZATION =================================================================================================
	
	# Calculating Ending Time
	timing_results['complete_analysis']['stop']=datetime.datetime.now()
	
	# Saving the timing results
	save_time_results(timing_results,optimization_input_parameters)
	
#=================================================================================================================================================================================
