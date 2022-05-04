#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the frequency analysis by sweeping frequency
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
	
	newpath =filename+'/Frequency_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Frequency_Analysis/Results/initial_circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in initial_circuit_parameters:
		f.write(str(param_name)+'\t'+str(initial_circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
def write_circuit_parameters(circuit_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	# Getting the filename
	
	newpath =filename+'/Frequency_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Frequency_Analysis/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each frequency iteration to a csv file
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,freq):
	filename=optimization_input_parameters['filename']['output']+'/Frequency_Analysis/Results/extracted_parameters.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(freq))
	for param_name in extracted_parameters:
		f.write(','+str(extracted_parameters[param_name]))	# Writing the values to the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Get the required data from the extracted parameters
def get_single_freq_extracted_parameters(extracted_parameters,freq):
	
	final_extracted_parameters={}
	for param in extracted_parameters:
		str_freq=str(freq)
		len_freq=len(str_freq)
		if str_freq==param[:len_freq]:
			final_extracted_parameters[param[len_freq+1:]]=extracted_parameters[param]
	
	return final_extracted_parameters


"""	
============================================================================================================================
------------------------------------- Output Functions ---------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the frequency analysis
def frequency_analysis(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['frequency_analysis']['run']=='NO':
		return

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Frequency Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['frequency_analysis']={}
	timing_results['frequency_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Frequency Analysis ***************************************************')
	
	cir.update_simulation_parameters(optimization_input_parameters['frequency_analysis']['simulation'])
	
	initial_circuit_parameters_initial=cir.get_initial_circuit_parameters()
	circuit_parameters_initial=cir.get_circuit_parameters()
	extracted_parameters_initial=cir.get_extracted_parameters()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Creating an array for frequency analysis
	freq_array=cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_list']
	
	# Writing the values to output files
	cir.run_circuit()
	write_initial_circuit_parameters(cir.initial_circuit_parameters,optimization_input_parameters)
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	
	# Storing the results
	extracted_parameters_iter={}
	
	# Performing the analysis
	for freq in freq_array:
		final_extracted_parameters=get_single_freq_extracted_parameters(cir.extracted_parameters,freq)
		update_extracted_parameters(final_extracted_parameters,optimization_input_parameters,freq)		# Writing the values to the output file
		extracted_parameters_iter[freq]=final_extracted_parameters.copy()

	# Restoring the value of initial extracted and circuit parameters
	cir.update_circuit_state(initial_circuit_parameters_initial,circuit_parameters_initial,extracted_parameters_initial)

	# Plotting the graphs
	file_directory=optimization_input_parameters['filename']['output']
	
	plot_type=optimization_input_parameters['frequency_analysis']['plot_type']
	plot_frequency_analysis(extracted_parameters_iter,file_directory,plot_type)

	# Storing the starting time
	timing_results['frequency_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Frequency Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()


"""
============================================================================================================================
------------------------------------- Plotting Functions -------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Frequency )
def plot_frequency_analysis(extracted_parameters_iter,file_directory,sweep_type):
	
	# Getting the file directory
	file_sub_directory=file_directory+'/Frequency_Analysis/Plots/'
	if not os.path.exists(file_sub_directory):
		os.makedirs(file_sub_directory)


	# Getting the data
	freq_array=[key for key in extracted_parameters_iter]
	parameters_list=[key for key in extracted_parameters_iter[freq_array[0]]]

	# Plotting the data
	for parameter in parameters_list:
		parameter_array=[extracted_parameters_iter[key][parameter] for key in freq_array]

		if sweep_type=='linear':
			figure()
			plot(freq_array,parameter_array)
			xlabel('Frequency')
			ylabel(parameter)
			grid()
			savefig(file_sub_directory+str(parameter)+'.pdf')
			close()
		
		else:
			figure()
			semilogx(freq_array,parameter_array)
			xlabel('Frequency')
			ylabel(parameter)
			grid()
			savefig(file_sub_directory+str(parameter)+'.pdf')
			close()


#===========================================================================================================================
