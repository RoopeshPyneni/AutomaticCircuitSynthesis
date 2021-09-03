#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: complete_optimization.py
File Description 	: This file will perform all the optimization steps by calling functions in different files

Functions structure in this file:
	--> complete_optimizations

COMPLETE
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
	fw.save_input_results(optimization_input_parameters)

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

