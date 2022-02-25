#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform all the optimization steps by calling functions in different files
"""

#===========================================================================================================================
import datetime
import optimization as op
import common_functions as cf # type: ignore
import CG_LNA.spectre as sp1
import CS_LNA.spectre as sp2
import os

# Analysis Files
import Analysis.sensitivity_analysis as sa
import Analysis.temperature_analysis as ta
import Analysis.process_analysis as pa
import Analysis.iip3_analysis as ia
import Analysis.frequency_analysis as fa
import Analysis.circuit_parameter_analysis as cpa


"""
===========================================================================================================================
------------------------------------ FILE WRITING FUNCTIONS ---------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that stores input data of the simulation ( Filenames, Output Conditions )
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
	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filenames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))

	# Saving Output Conditions
	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output Conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	for name in optimization_input_parameters['output_conditions']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['output_conditions'][name],3))

	f.close()

#-----------------------------------------------------------------
# Function that creates a file that stores output data of the simulation
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
	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filenames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))
	f.write('\n\n\n')

	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the MOS File Calculations
def save_mos_results(mos_parameters,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+'/output_data.txt'
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOS Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

	for param_name in mos_parameters:
		f.write('\n'+str(param_name)+': '+cf.num_trunc(mos_parameters[param_name],3))
	
	f.close()

#-----------------------------------------------------------------
# Function that stores the time results
def save_time_results(timing_results,optimization_input_parameters):
	
	# Opening the filename
	filename=optimization_input_parameters['filename']['output']	
	filename=filename+str('/output_data.txt')
	f=open(filename,'a')

	# Storing the timing results
	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Timing Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	
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


"""
===========================================================================================================================
------------------------------------ MAIN PROGRAM -------------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that performs the complete optimization process
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

	#=============================== MOSFET PARAMETERS =============================================

	# Writing the MOSFET File Location to .scs file
	if name=='CG_LNA':
		cir=sp1.Circuit(circuit_initialization_parameters)
	else:
		cir=sp2.Circuit(circuit_initialization_parameters)
	save_mos_results(cir.mos_parameters,optimization_input_parameters)

	#=============================== PRE OPTIMIZATION ==============================================

	cir.pre_optimization(optimization_input_parameters,timing_results)

	#=============================== OPTIMIZATION ==================================================

	op.main_opt(cir,optimization_input_parameters,timing_results)

	#=============================== SENSITIVITY ANALYSIS ==========================================

	sa.sensitivity_analysis(cir,optimization_input_parameters,timing_results)

	#=============================== TEMPERATURE ANALYSIS ==========================================

	ta.temperature_analysis(cir,optimization_input_parameters,timing_results)

	#=============================== PROCESS ANALYSIS ==============================================

	pa.process_analysis(cir,optimization_input_parameters,timing_results)

	#=============================== IIP3 ANALYSIS =================================================

	ia.iip3_analysis(cir,optimization_input_parameters,timing_results)

	#=============================== FREQUENCY ANALYSIS ============================================

	fa.frequency_analysis(cir,optimization_input_parameters,timing_results)
	
	#=============================== CIRCUIT PARAMETER ANALYSIS ============================================

	cpa.circuit_parameter_analysis(cir,optimization_input_parameters,timing_results)
	
	#=============================== AFTER OPTIMIZATION ============================================
	
	# Calculating Ending Time
	timing_results['complete_analysis']['stop']=datetime.datetime.now()
	
	# Saving the timing results
	save_time_results(timing_results,optimization_input_parameters)


#=================================================================================================================================================================================
