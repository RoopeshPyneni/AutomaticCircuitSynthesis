#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the sensitivity analysis by varying each circuit parameter by 1%
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf # type: ignore
from matplotlib import pylab
from pylab import *
#===========================================================================================================================


"""
============================================================================================================================
------------------------------------- Defining the functions ---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
def write_initial_circuit_parameters(initial_circuit_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	# Getting the filename
	
	newpath =filename+'/Sensitivity_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Sensitivity_Analysis/Results/initial_circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in initial_circuit_parameters:
		f.write(str(param_name)+'\t'+str(initial_circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
def write_circuit_parameters(circuit_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	# Getting the filename
	
	newpath =filename+'/Sensitivity_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Sensitivity_Analysis/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/extracted_parameters_sensitivity.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Circuit Parameter')
	for param_name in extracted_parameters:
		f.write(','+param_name)	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
def update_extracted_parameters(param,extracted_parameters_sensitivity,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/extracted_parameters_sensitivity.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(param))
	for param_name in extracted_parameters_sensitivity:
		f.write(','+str(extracted_parameters_sensitivity[param_name]))	# Writing the values to the file
	f.write('\n')
	f.close()
	


"""
============================================================================================================================
------------------------------------- Output Functions ---------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the temperature analysis
def sensitivity_analysis(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['sensitivity_analysis']['run']=='NO':
		return

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Sensitivity Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['sensitivity_analysis']={}
	timing_results['sensitivity_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Sensitivity Analysis ***************************************************')
	
	cir.update_simulation_parameters(optimization_input_parameters['sensitivity_analysis']['simulation'])
	
	initial_circuit_parameters_initial=cir.get_initial_circuit_parameters()
	circuit_parameters_initial=cir.get_circuit_parameters()
	extracted_parameters_initial=cir.get_extracted_parameters()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Writing the values to output files
	cir.run_circuit()
	write_initial_circuit_parameters(cir.initial_circuit_parameters,optimization_input_parameters)
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters)

	extracted_parameters_base=cir.extracted_parameters.copy()
	
	# Performing the analysis
	for param in initial_circuit_parameters_initial:
		initial_circuit_parameters_current=initial_circuit_parameters_initial.copy()
		initial_circuit_parameters_current[param]*=1.01
		cir.update_circuit(initial_circuit_parameters_current)
		extracted_parameters_sensitivity={}
		for ext_param in extracted_parameters_base:
			if extracted_parameters_base[ext_param]==0:
				extracted_parameters_sensitivity[ext_param]=0
			else:
				extracted_parameters_sensitivity[ext_param]=100*(cir.extracted_parameters[ext_param]-extracted_parameters_base[ext_param])/(extracted_parameters_base[ext_param])
		update_extracted_parameters(param,extracted_parameters_sensitivity,optimization_input_parameters)		# Writing the values to the output file
		
	# Restoring the value of initial extracted and circuit parameters
	cir.update_circuit_state(initial_circuit_parameters_initial,circuit_parameters_initial,extracted_parameters_initial)

	# Storing the starting time
	timing_results['sensitivity_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Sensitivity Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	
	
#===========================================================================================================================

# Plotting the graphs
#file_directory=optimization_input_parameters['filename']['output']
#spec_current=cir.circuit_parameters['Io']
#plot_sensitivity_analysis(extracted_parameters_iter,file_directory,spec_current)
