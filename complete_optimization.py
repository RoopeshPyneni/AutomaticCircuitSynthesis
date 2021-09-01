#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Full Optimization Code:
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
#------------------------------------ MOSFET EXTRACTION --------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that extracts the MOSFET File Parameeters
def calculate_mos_parameters(optimization_input_parameters):
	
	# Setting Lmin and Vdd
	Lmin=optimization_input_parameters['MOS']['Lmin']
	vdd=optimization_input_parameters['MOS']['Vdd']

	# Extracting From File
	mos_file_parameters = {'un':0,'cox':0,'vt':0,'Lmin':Lmin,'vdd':vdd}
	mos_file_parameters=sp.extract_mosfet_param(optimization_input_parameters,mos_file_parameters)
	mos_parameters=mos_file_parameters.copy()

	# Printing the MOSFET Parameters
	cf.print_MOS_parameters(mos_parameters)

	# Storing the results
	fw.save_mos_results(mos_parameters,optimization_input_parameters)

	return mos_parameters




#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------

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
	
	print('************************************************************************************************************')
	print('*********************************** MOSFET Parameters ******************************************************')

	# Extracting the MOSFET Parameters from the MOS file
	mos_parameters=calculate_mos_parameters(optimization_input_parameters)

	# Writing the MOSFET File Location to .scs file
	sp.write_MOS_parameters(optimization_input_parameters)

	#======================================================== PRE OPTIMIZATION ===================================================================================================

	print('************************************************************************************************************')
	print('*********************************** Pre Optimization *******************************************************')
	
	# Running pre optimization
	circuit_parameters,extracted_parameters=pr.pre_optimization(mos_parameters,optimization_input_parameters,timing_results)

	#======================================================== OPTIMIZATION =======================================================================================================

	print('************************************************************************************************************')
	print('*********************************** Main Optimization ******************************************************')

	# Running the main optimization
	circuit_parameters,extracted_parameters=op.main_opt(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)

	#======================================================== TEMPERATURE ANALYSIS ===============================================================================================

	print('************************************************************************************************************')
	print('*********************************** Temperature Analysis ***************************************************')

	# Running the temperature analysis
	circuit_parameters,extracted_parameters=ta.temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)

	#======================================================== TEMPERATURE ANALYSIS 2 ==============================================================================================

	print('************************************************************************************************************')
	print('*********************************** Temperature Analysis ***************************************************')

	# Running the temperature analysis
	circuit_parameters,extracted_parameters=ta2.temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results)
	
	#======================================================== AFTER OPTIMIZATION =================================================================================================
	
	# Calculating Ending Time
	timing_results['complete_analysis']['stop']=datetime.datetime.now()
	
	fw.save_time_results(timing_results,optimization_input_parameters)
	
	#=============================================================================================================================================================================

