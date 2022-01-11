#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the frequency analysis by sweeping frequency
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf
from matplotlib import pylab
from pylab import *
#===========================================================================================================================


#===========================================================================================================================
#------------------------------------Defining the functions -----------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
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
# Writing the header row for extracted parameters to a csv file
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Frequency_Analysis/Results/extracted_parameters.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Frequency')
	for param_name in extracted_parameters:
		f.write(','+param_name)	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each frequency iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,freq):
	filename=optimization_input_parameters['filename']['output']+'/Frequency_Analysis/Results/extracted_parameters.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(freq))
	for param_name in extracted_parameters:
		f.write(','+str(extracted_parameters[param_name]))	# Writing the values to the file
	f.write('\n')
	f.close()
	


		
#===========================================================================================================================
#--------------------------------------------Output Functions---------------------------------------------------------------
	
#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the frequency analysis
# Input: circuit_parameters, extracted_parameters, optimization_input_parameters
# Output: circuit_parameters, extracted_parameters
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
	
	initial_extracted_parameters=cir.extracted_parameters.copy()
	initial_circuit_parameters=cir.circuit_parameters.copy()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Creating an array for frequency analysis
	start_freq=optimization_input_parameters['frequency_analysis']['start_freq']
	stop_freq=optimization_input_parameters['frequency_analysis']['stop_freq']
	n_freq=optimization_input_parameters['frequency_analysis']['n_freq']
	sweep_type=optimization_input_parameters['frequency_analysis']['sweep_type']
	if sweep_type=='linear':
		freq_array=np.linspace(start_freq,stop_freq,n_freq)
	else:
		freq_array=np.logspace(np.log10(start_freq),np.log10(stop_freq),n_freq)
	
	# Writing the values to output files
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters)

	# Getting the initial value of frequency
	initial_frequency=cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']

	# Storing the results
	extracted_parameters_iter={}
	
	# Performing the analysis
	for freq in freq_array:
		cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=freq
		cir.run_circuit()
		update_extracted_parameters(cir.extracted_parameters,optimization_input_parameters,freq)		# Writing the values to the output file
		extracted_parameters_iter[freq]=cir.extracted_parameters.copy()

	# Restoring the value of initial extracted and circuit parameters
	cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=initial_frequency
	cir.run_circuit()

	# Plotting the graphs
	file_directory=optimization_input_parameters['filename']['output']
	spec_current=cir.circuit_parameters['Io']
	plot_frequency_analysis(extracted_parameters_iter,file_directory,sweep_type)

	# Storing the starting time
	timing_results['frequency_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Frequency Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	
#===========================================================================================================================


"""
====================================================================================================================================================================
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Frequency )
def plot_frequency_analysis(extracted_parameters_iter,file_directory,sweep_type):
	
	# Getting the file directory
	file_sub_directory=file_directory+'/Frequency_Analysis/Plots/'

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
