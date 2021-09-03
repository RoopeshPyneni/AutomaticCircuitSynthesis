#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: complete_optimization.py
File Description 	: This file will perform all the optimization steps by calling functions in different files

Functions structure in this file:
	--> complete_optimization

"""

#===========================================================================================================================
import datetime
import optimization as op
import common_functions as cf
import pre_optimization as pr
import temperature_analysis as ta
import temperature_analysis_2 as ta2
import spectre as sp
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

	# Saving MOS Parameters
	f.write('\n\n---------------------- MOS Parameters -----------------------')
	f.write('\nMOSFET File	:'+str(optimization_input_parameters['MOS']['filename']))
	f.write('\nMOS Type 	:'+str(optimization_input_parameters['MOS']['Type']))
	f.write('\nVdd      	:'+str(optimization_input_parameters['MOS']['Vdd']))
	f.write('\nLmin     	:'+str(optimization_input_parameters['MOS']['Lmin']))

	# Saving Output Conditions
	f.write('\n\n---------------------- Output Conditions -----------------------')
	for name in optimization_input_parameters['output_conditions']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['output_conditions'][name],3))

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
	
	# Circuit Writing List for simulation
	for name in optimization_input_parameters['simulation']['cir_writing_dict']:
		f.write('\n'+str(name)+': '+str(optimization_input_parameters['simulation']['cir_writing_dict'][name]))

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

#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that performs the complete optimization process
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def complete_optimization(optimization_input_parameters):

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
	sp.write_MOS_parameters(optimization_input_parameters)

	#======================================================== PRE OPTIMIZATION ===================================================================================================

	circuit_parameters,extracted_parameters=pr.pre_optimization(optimization_input_parameters,timing_results)

	#======================================================== OPTIMIZATION =======================================================================================================

	circuit_parameters,extracted_parameters=op.main_opt(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)

	#======================================================== TEMPERATURE ANALYSIS ===============================================================================================

	circuit_parameters,extracted_parameters=ta.temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)

	#======================================================== TEMPERATURE ANALYSIS 2 =============================================================================================

	circuit_parameters,extracted_parameters=ta2.temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)
	
	#======================================================== AFTER OPTIMIZATION =================================================================================================
	
	# Calculating Ending Time
	timing_results['complete_analysis']['stop']=datetime.datetime.now()
	
	# Saving the timing results
	save_time_results(timing_results,optimization_input_parameters)
	
#=================================================================================================================================================================================