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
import file_write as fw
import optimization as op
import common_functions as cf
import pre_optimization as pr
import temperature_analysis as ta
import temperature_analysis_2 as ta2
import spectre as sp
import os

#===========================================================================================================================
#------------------------------------ File Storage Code --------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that stores input data of the simulation
def save_input_results_initial(optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	filename=filename+str('/input_data.txt')
	f=open(filename,'w')

	# Saving Output Conditions
	f.write('\n\n---------------------- Output Conditions -----------------------')
	for name in optimization_input_parameters['output_conditions']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['output_conditions'][name],3))

	# Saving MOS Parameters
	f.write('\n\n---------------------- MOS Parameters -----------------------')
	f.write('\nMOSFET File	:'+str(optimization_input_parameters['MOS']['filename']))
	f.write('\nMOS Type 	:'+str(optimization_input_parameters['MOS']['Type']))
	f.write('\nVdd      	:'+str(optimization_input_parameters['MOS']['Vdd']))
	f.write('\nLmin     	:'+str(optimization_input_parameters['MOS']['Lmin']))

	# Saving Filenames
	f.write('\n\n---------------------- Filenames -----------------------')
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))

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

	# Saving the optimization input results
	save_input_results_initial(optimization_input_parameters)

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

	#======================================================== TEMPERATURE ANALYSIS 2 ==============================================================================================

	circuit_parameters,extracted_parameters=ta2.temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)
	
	#======================================================== AFTER OPTIMIZATION =================================================================================================
	
	# Calculating Ending Time
	timing_results['complete_analysis']['stop']=datetime.datetime.now()
	
	fw.save_time_results(timing_results,optimization_input_parameters)
	
	#=============================================================================================================================================================================

