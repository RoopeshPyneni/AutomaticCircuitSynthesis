#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Temperature Analysis:
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
	
	newpath =filename+'/Temperature_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Temperature_Analysis/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis/Results/extracted_parameters.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Temp,')
	for param_name in extracted_parameters:
		f.write(param_name+',')	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(current_temp,extracted_parameters,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis/Results/extracted_parameters.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(current_temp)+',')
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
def temperature_analysis(circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['temperature_analysis']['run']=='NO':
		return circuit_parameters,extracted_parameters

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['temperature_analysis']={}
	timing_results['temperature_analysis']['start']=datetime.datetime.now()

	initial_extracted_parameters=extracted_parameters.copy()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Creating an array for temperature analysis
	start_temp=optimization_input_parameters['temperature_analysis']['start_temp']
	stop_temp=optimization_input_parameters['temperature_analysis']['stop_temp']
	n_temp=optimization_input_parameters['temperature_analysis']['n_temp']
	temp_array=np.linspace(start_temp,stop_temp,n_temp)

	# Writing the values to output files
	write_circuit_parameters(circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters)

	# Performing the analysis
	for current_temp in temp_array:
		optimization_input_parameters['simulation']['parameters_list']['cir_temp']=current_temp			# Writing the temperature value to the netlist file
		extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)			# Extracting the parameters
		update_extracted_parameters(current_temp,extracted_parameters,optimization_input_parameters)	# Writing the values to the output file
		extracted_parameters_iter[current_temp]=extracted_parameters.copy()

	extracted_parameters=initial_extracted_parameters.copy()
	optimization_input_parameters['simulation']['parameters_list']['cir_temp']=optimization_input_parameters['simulation']['std_temp']	# Writing the temperature value to the netlist file

	# Plotting the graphs
	file_directory=optimization_input_parameters['filename']['output']
	plot_temp_analysis(extracted_parameters_iter,file_directory)

	# Storing the ending time
	timing_results['temperature_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	
	return circuit_parameters,extracted_parameters
	
#===========================================================================================================================

"""
================================================================================================================================================================================================================
"""


#-----------------------------------------------------------------------------------------------
# Function to extract the data from temperature analysis csv file
def extract_temp_analysis(extracted_parameters_iter):
	
	# Assigning the array to store the loss variables
	param_name_list=[]
	
	# Creating variables for the counting no of lines and variables
	n_iter=0
	flag=0
	for temp in extracted_parameters_iter:
		n_iter+=1
		if flag==0:
			temp_initial=temp
			flag=1
	
	n_param=1
	param_name_list.append('Temp')
	for param_name in extracted_parameters_iter[temp_initial]:
		n_param+=1
		param_name_list.append(param_name)
	
	# Creating array to store the values 
	extracted_array=np.zeros((n_iter,n_param),dtype=float)
	
	# Extracting the values of the variables
	i=0
	for temp in extracted_parameters_iter:
		j=1
		# Storing the variables
		extracted_array[i,0]=temp
		for param_name in extracted_parameters_iter[temp]:
			extracted_array[i,j]=extracted_parameters_iter[temp][param_name]
			j+=1
		i+=1
		
	return extracted_array,param_name_list

#-----------------------------------------------------------------------------------------------
# Plotting results from temperature analysis
def plot_temp_analysis(extracted_parameters_iter,file_directory):
	
	extracted_array,param_name=extract_temp_analysis(extracted_parameters_iter)

	pathname=file_directory+'/Temperature_Analysis/Plots/'
	
	if not os.path.exists(pathname):
		os.makedirs(pathname)
    		
	n_param=len(param_name)

	arrX=extracted_array[:,0]

	output_conditions_dict={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}

	# Figure 1
	for i in range(n_param-1):
		figure()
		plot(arrX,extracted_array[:,i+1],'g',label=param_name[i+1])
		if param_name[i+1] in output_conditions_dict:
			arrY=np.ones(len(arrX),dtype=float)
			arrY=arrY*output_conditions_dict[param_name[i+1]]
			plot(arrX,arrY,'r',label='Output Condition')
		xlabel('Temperature')
		ylabel(param_name[i+1])
		legend()
		grid()
		savefig(pathname+param_name[i+1]+'.pdf')
		close()

#===========================================================================================================================




		
	


