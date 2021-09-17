#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: process_analysis.py
File Description 	: This file will perform the process analysis by changing the process and measuring the values

Functions structure in this file:
	--> process_analysis
		--> write_circuit_parameters
		--> write_extracted_parameters_initial
		--> update_extracted_parameters
		
		--> plot_process_analysis
			--> extract_process_analysis

COMPLETE
"""

#===========================================================================================================================
import numpy as np
import spectre as sp
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
	
	newpath =filename+'/Process_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Process_Analysis/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Process_Analysis/Results/extracted_parameters.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Process,')
	for param_name in extracted_parameters:
		f.write(param_name+',')	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each process iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,process_corner):
	filename=optimization_input_parameters['filename']['output']+'/Process_Analysis/Results/extracted_parameters.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(process_corner)+',')
	for param_name in extracted_parameters:
		f.write(str(extracted_parameters[param_name])+',')	# Writing the values to the file
	f.write('\n')
	f.close()
	


		
#===========================================================================================================================
#--------------------------------------------Output Functions---------------------------------------------------------------
	
#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the temperature analysis
# Input: circuit_parameters, extracted_parameters, optimization_input_parameters
# Output: circuit_parameters, extracted_parameters
def process_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['process_analysis']['run']=='NO':
		return circuit_parameters,extracted_parameters

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Process Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['process_analysis']={}
	timing_results['process_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('************************************* Process Analysis *****************************************************')
	
	cf.write_simulation_parameters(optimization_input_parameters,'process_analysis',0)

	initial_extracted_parameters=extracted_parameters.copy()
	initial_circuit_parameters=circuit_parameters.copy()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Writing the values to output files
	write_circuit_parameters(circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters)

	# Creating the list of process corners
	process_corner_list=['tt','ff','ss']
	
	# Performing the analysis
	for process_corner in process_corner_list:
		extracted_parameters_iter[process_corner]={}
		optimization_input_parameters['simulation']['process_corner']=process_corner
		sp.write_simulation_parameters(optimization_input_parameters)
		extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)				# Extracting the parameters
		update_extracted_parameters(extracted_parameters,optimization_input_parameters,process_corner)		# Writing the values to the output file
		
	# Restoring the value of initial extracted and circuit parameters
	extracted_parameters=initial_extracted_parameters.copy()
	circuit_parameters=initial_circuit_parameters.copy()

	# Plotting the graphs
	#file_directory=optimization_input_parameters['filename']['output']
	#plot_process_analysis(extracted_parameters_iter,file_directory)

	# Storing the stoppingg time
	timing_results['process_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Process Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	
	return circuit_parameters,extracted_parameters
	
#===========================================================================================================================


"""
====================================================================================================================================================================
"""

"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
# Inputs  : extracted_parameters_iter
# Outputs : extracted_matrix, temp_array, current_array, param_array
def plot_process_analysis(extracted_parameters_iter,file_directory,spec_current):
	
	# Creating a folder to store the results
	pathname=file_directory+'/Process_Analysis/Plots/'
	if not os.path.exists(pathname):
		os.makedirs(pathname)
	
	# Creating output conditions dictionary
	output_conditions={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}
	#colour_dict={'iip3_dbm':'r','nf_db':'g','s11_db':'b','gain_db':'m'}
	colour_dict={0:'green',1:'blue',2:'red',3:'cyan',4:'lime',5:'violet'}

	extracted_matrix,process_array,param_array=extract_process_analysis(extracted_parameters_iter)

	# Getting the lengths of the arrays
	n_process=len(process_array)
	n_param=len(param_array)

#-----------------------------------------------------------------------------------------------
# Function to extract the data from extracted_parameters_iter dictionary and store it in the form of a matrix
# Inputs  : extracted_parameters_iter
# Outputs : extracted_matrix, temp_array, current_array, param_array
def extract_process_analysis(extracted_parameters_iter):
	
	# Assigning the array to store the loss variables
	process_array=[]
	param_array=[]
	
	# Creating variables for the counting no of lines and variables
	n_process=0
	n_param=0
	
	# Storing the value of temp, current, and param arrays
	for process in extracted_parameters_iter:
		process_array.append(process)
		n_process+=1
		final_process=process

	for param in extracted_parameters_iter[final_process]:
		param_array.append(param)
		n_param+=1
		
	# Creating array to store the values 
	extracted_matrix=np.zeros((n_process,n_param),dtype=float)
	
	# Extracting the values of the variables
	i=0
	for process in extracted_parameters_iter:
		j=0
		for param in extracted_parameters_iter[process]:
			extracted_matrix[i,j]=extracted_parameters_iter[process][param]
			j+=1
		i+=1
	
	return extracted_matrix,process_array,param_array

"""

#===========================================================================================================================
