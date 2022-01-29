#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the circuit parameter sweep analysis by sweeping a circuit parameter
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
def write_circuit_parameters(circuit_parameters,optimization_input_parameters,i):
	
	filename=optimization_input_parameters['filename']['output']	# Getting the filename
	
	newpath =filename+'/Circuit_Parameter_Analysis/'+str(i)+'/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Circuit_Parameter_Analysis/'+str(i)+'/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters,parameter_name,i):
	
	filename=optimization_input_parameters['filename']['output']+'/Circuit_Parameter_Analysis/'+str(i)+'/Results/extracted_parameters.csv'	# Getting the filename

	f=open(filename,'w')
	f.write(str(parameter_name))
	for param_name in extracted_parameters:
		f.write(','+param_name)	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each iteration to a csv file
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,parameter_value,i):
	filename=optimization_input_parameters['filename']['output']+'/Circuit_Parameter_Analysis/'+str(i)+'/Results/extracted_parameters.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(parameter_value))
	for param_name in extracted_parameters:
		f.write(','+str(extracted_parameters[param_name]))	# Writing the values to the file
	f.write('\n')
	f.close()


"""	
============================================================================================================================
------------------------------------- Output Functions ---------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the circuit parameter analysis
def circuit_parameter_analysis(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['circuit_parameter_analysis']['run']=='NO':
		return
	n_runs=optimization_input_parameters['circuit_parameter_analysis']['n_runs']

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Circuit Parameter Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['circuit_parameter_analysis']={}
	timing_results['circuit_parameter_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Circuit Parameter Analysis *********************************************')

	# Updating the simulation parameters
	cir.update_simulation_parameters(optimization_input_parameters['circuit_parameter_analysis']['simulation'])

	for i in range(n_runs):	
		
		# Getting the initial value of parameters
		initial_extracted_parameters=cir.extracted_parameters.copy()
		initial_circuit_parameters=cir.circuit_parameters.copy()
		
		# Creating an array for circuit parameter analysis
		parameter_name=optimization_input_parameters['circuit_parameter_analysis'][i]['parameter_name']
		parameter_select_type=optimization_input_parameters['circuit_parameter_analysis'][i]['parameter_select_type']
		start=optimization_input_parameters['circuit_parameter_analysis'][i]['start']
		stop=optimization_input_parameters['circuit_parameter_analysis'][i]['stop']
		n_value=optimization_input_parameters['circuit_parameter_analysis'][i]['n_value']
		sweep_type=optimization_input_parameters['circuit_parameter_analysis'][i]['sweep_type']

		if parameter_select_type=='relative':
			start*=cir.circuit_parameters[parameter_name]
			stop*=cir.circuit_parameters[parameter_name]

		if sweep_type=='linear':
			parameter_array=np.linspace(start,stop,n_value)
		else:
			parameter_array=np.logspace(np.log10(start),np.log10(stop),n_value)
		
		# Writing the values to output files
		write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters,i)
		write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters,parameter_name,i)

		# Storing the results
		extracted_parameters_iter={}
		
		# Performing the analysis
		for value in parameter_array:
			cir.circuit_parameters[parameter_name]=value
			cir.run_circuit()
			update_extracted_parameters(cir.extracted_parameters,optimization_input_parameters,value,i)		# Writing the values to the output file
			extracted_parameters_iter[value]=cir.extracted_parameters.copy()

		# Restoring the value of initial extracted and circuit parameters
		cir.circuit_parameters=initial_circuit_parameters.copy()
		cir.extracted_parameters=initial_extracted_parameters.copy()

		# Plotting the graphs
		file_directory=optimization_input_parameters['filename']['output']
		plot_circuit_parameter_analysis(extracted_parameters_iter,file_directory,sweep_type,parameter_name,i)

	# Storing the ending time
	timing_results['circuit_parameter_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Circuit Parameter Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()


"""
============================================================================================================================
------------------------------------- Plotting Functions -------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Output Parameters vs Input Parameter )
def plot_circuit_parameter_analysis(extracted_parameters_iter,file_directory,sweep_type,parameter_name,i):
	
	# Getting the file directory
	file_sub_directory=file_directory+'/Circuit_Parameter_Analysis/'+str(i)+'/Plots/'
	if not os.path.exists(file_sub_directory):
		os.makedirs(file_sub_directory)

	# Getting the data
	circuit_parameter_array=[key for key in extracted_parameters_iter]
	parameters_list=[key for key in extracted_parameters_iter[circuit_parameter_array[0]]]

	# Plotting the data
	for parameter in parameters_list:
		parameter_array=[extracted_parameters_iter[key][parameter] for key in circuit_parameter_array]

		if sweep_type=='linear':
			figure()
			plot(circuit_parameter_array,parameter_array)
			xlabel(parameter_name)
			ylabel(parameter)
			grid()
			savefig(file_sub_directory+str(parameter)+'.pdf')
			close()
		
		else:
			figure()
			semilogx(circuit_parameter_array,parameter_array)
			xlabel(parameter_name)
			ylabel(parameter)
			grid()
			savefig(file_sub_directory+str(parameter)+'.pdf')
			close()


#===========================================================================================================================
