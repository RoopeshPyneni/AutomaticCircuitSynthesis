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
	
	newpath =filename+'/Temperature_Analysis_2/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Temperature_Analysis_2/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters,temp):
	
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis_2/Results/extracted_parameters_'+str(temp)+'.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Current Mirror Io,')
	for param_name in extracted_parameters:
		f.write(param_name+',')	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,temp,current):
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis_2/Results/extracted_parameters_'+str(temp)+'.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(current)+',')
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
	
	if optimization_input_parameters['temperature_analysis_2']['run']=='NO':
		return circuit_parameters,extracted_parameters

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis 2 Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['temperature_analysis_2']={}
	timing_results['temperature_analysis_2']['start']=datetime.datetime.now()

	initial_extracted_parameters=extracted_parameters.copy()
	initial_circuit_parameters=circuit_parameters.copy()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Creating an array for temperature analysis
	start_temp=optimization_input_parameters['temperature_analysis_2']['start_temp']
	stop_temp=optimization_input_parameters['temperature_analysis_2']['stop_temp']
	n_temp=optimization_input_parameters['temperature_analysis_2']['n_temp']
	temp_array=np.linspace(start_temp,stop_temp,n_temp)
	
	# Creating an array for Io sweep
	room_temp_current_log=np.log10(circuit_parameters['Io'])
	start_current=np.log10(optimization_input_parameters['temperature_analysis_2']['start_current'])
	stop_current=np.log10(optimization_input_parameters['temperature_analysis_2']['stop_current'])
	n_current=optimization_input_parameters['temperature_analysis_2']['n_current']
	current_array=np.logspace(room_temp_current_log+start_current,room_temp_current_log+stop_current,n_current)
	
	# Writing the values to output files
	write_circuit_parameters(circuit_parameters,optimization_input_parameters)
	
	# Performing the analysis
	for temp in temp_array:
		extracted_parameters_iter[temp]={}
		write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters,temp)
		for current in current_array:
			circuit_parameters['Io']=current
			optimization_input_parameters['simulation']['parameters_list']['cir_temp']=temp				# Writing the temperature value to the netlist file
			extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)			# Extracting the parameters
			update_extracted_parameters(extracted_parameters,optimization_input_parameters,temp,current)		# Writing the values to the output file
			extracted_parameters_iter[temp][current]=extracted_parameters.copy()

	# Restoring the value of initial extracted and circuit parameters
	extracted_parameters=initial_extracted_parameters.copy()
	circuit_parameters=initial_circuit_parameters.copy()

	optimization_input_parameters['simulation']['parameters_list']['cir_temp']=optimization_input_parameters['simulation']['std_temp']	# Writing the temperature value to the netlist file

	# Plotting the graphs
	file_directory=optimization_input_parameters['filename']['output']
	#plot_temp_analysis(extracted_parameters_iter,file_directory)
	
	spec_current=circuit_parameters['Io']
	plot_temp_analysis_2(extracted_parameters_iter,file_directory,spec_current)

	# Storing the starting time
	timing_results['temperature_analysis_2']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis 2 End\n Time : '+str(datetime.datetime.now())+'\n\n')
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
	for current in extracted_parameters_iter:
		n_iter+=1
		if flag==0:
			current_initial=current
			flag=1
	
	n_param=1
	param_name_list.append('Io')
	for param_name in extracted_parameters_iter[current_initial]:
		n_param+=1
		param_name_list.append(param_name)
	
	# Creating array to store the values 
	extracted_array=np.zeros((n_iter,n_param),dtype=float)
	
	# Extracting the values of the variables
	i=0
	for current in extracted_parameters_iter:
		j=1
		# Storing the variables
		extracted_array[i,0]=current
		for param_name in extracted_parameters_iter[current]:
			extracted_array[i,j]=extracted_parameters_iter[current][param_name]
			j+=1
		i+=1
		
	return extracted_array,param_name_list

#-----------------------------------------------------------------------------------------------
# Plotting results from temperature analysis
def plot_single_temp(extracted_parameters_iter,file_directory):
	
	extracted_array,param_name=extract_temp_analysis(extracted_parameters_iter)
	
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
    		
	n_param=len(param_name)

	arrX=extracted_array[:,0]

	output_conditions_dict={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}

	# Figure 1
	for i in range(n_param-1):
		figure()
		semilogx(arrX,extracted_array[:,i+1],'g',label=param_name[i+1])
		if param_name[i+1] in output_conditions_dict:
			arrY=np.ones(len(arrX),dtype=float)
			arrY=arrY*output_conditions_dict[param_name[i+1]]
			semilogx(arrX,arrY,'r',label='Output Condition')
		xlabel('Io')
		ylabel(param_name[i+1])
		legend()
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(file_directory+param_name[i+1]+'.pdf')
		close()
	
	# Figure 2
	colour_dict={'iip3_dbm':'r','nf_db':'g','s11_db':'b','gain_db':'m'}
	figure()
	for i in range(n_param-1):
		if param_name[i+1] in output_conditions_dict:
			semilogx(arrX,extracted_array[:,i+1],colour_dict[param_name[i+1]],linestyle='-',label=param_name[i+1])
			arrY=np.ones(len(arrX),dtype=float)
			arrY=arrY*output_conditions_dict[param_name[i+1]]
			semilogx(arrX,arrY,colour_dict[param_name[i+1]],linestyle='--',label=param_name[i+1]+' Output Condition')
	xlabel('Io')
	ylabel('Output Parameters')
	legend()
	grid(b=True,which='major',color='#666666')
	grid(b=True,which='minor',color='#999999')
	savefig(file_directory+'all.pdf')
	close()

#-----------------------------------------------------------------------------------------------
# Plotting results from temperature analysis
def plot_temp_analysis(extracted_parameters_iter,file_directory):
	for temp in extracted_parameters_iter:
		file_sub_directory=file_directory+'/Temperature_Analysis_2/Plots/'+str(temp)+'/'
		plot_single_temp(extracted_parameters_iter[temp],file_sub_directory)

"""
====================================================================================================================================================================
"""

#-----------------------------------------------------------------------------------------------
# Function to extract the data from extracted_parameters_iter dictionary and store it in the form of a matrix
def extract_temp_analysis_2(extracted_parameters_iter):
	
	# Assigning the array to store the loss variables
	temp_array=[]
	current_array=[]
	param_array=[]
	
	# Creating variables for the counting no of lines and variables
	n_temp=0
	n_current=0
	n_param=0
	
	# Storing the value of temp, current, and param arrays
	for temp in extracted_parameters_iter:
		temp_array.append(temp)
		n_temp+=1
		final_temp=temp

	for current in extracted_parameters_iter[final_temp]:
		current_array.append(current)
		n_current+=1
		final_current=current

	for param in extracted_parameters_iter[final_temp][final_current]:
		param_array.append(param)
		n_param+=1
		
	# Creating array to store the values 
	extracted_matrix=np.zeros((n_temp,n_current,n_param),dtype=float)
	
	# Extracting the values of the variables
	i=0
	for temp in extracted_parameters_iter:
		j=0
		for current in extracted_parameters_iter[temp]:
			k=0
			for param in extracted_parameters_iter[temp][current]:
				extracted_matrix[i,j,k]=extracted_parameters_iter[temp][current][param]
				k+=1
			j+=1
		i+=1
	
	return extracted_matrix,temp_array,current_array,param_array

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
def plot_temp_analysis_2(extracted_parameters_iter,file_directory,spec_current):
	
	file_sub_directory=file_directory+'/Temperature_Analysis_2/Plots/'
	extracted_matrix,temp_array,current_array,param_array=extract_temp_analysis_2(extracted_parameters_iter)

	output_conditions={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}
	colour_dict={'iip3_dbm':'r','nf_db':'g','s11_db':'b','gain_db':'m'}


	plot_param_vs_current(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_current/',output_conditions,colour_dict,spec_current)
	plot_param_vs_temperature(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_temperature/',output_conditions,colour_dict,spec_current)

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs current at different temperatures )
def plot_param_vs_current(extracted_matrix,temp_array,current_array,param_array,file_sub_directory,output_conditions,colour_dict,spec_current):
	
	# Finding the length of each array
	n_temp=len(temp_array)
	n_current=len(current_array)
	n_param=len(param_array)

	# First, we will plot Parameters vs Current and each plot will contain a single temperature
	for i in range(n_temp):
		# Creating the directory
		pathname=file_sub_directory+str(temp_array[i])+'/'
		if not os.path.exists(pathname):
			os.makedirs(pathname)

		for k in range(n_param):
			
			# Plot for single parameter
			figure()
			
			# Plotting the main parameter
			semilogx(current_array,extracted_matrix[i,:,k],color='g',label=param_array[k])
			
			# Calculating the maximum and minimum value of the output
			max_val=np.max(extracted_matrix[i,:,k])
			min_val=np.min(extracted_matrix[i,:,k])
			if param_array[k] in output_conditions:
				arrY=output_conditions[param_array[k]]*np.ones(n_current,dtype=float)
				semilogx(current_array,arrY,color='r',label='Output Spec')
				if output_conditions[param_array[k]]>max_val:
					max_val=output_conditions[param_array[k]]
				if output_conditions[param_array[k]]<min_val:
					min_val=output_conditions[param_array[k]]

			# Plotting a vertical line at spec_current
			arrX=spec_current*np.ones(100,dtype=float)
			arrY=np.linspace(min_val,max_val,100)
			semilogx(arrX,arrY,color='black',linestyle='--')

			# Other Plotting Parameters
			xlabel('Io')
			ylabel(param_array[k])
			legend()
			grid(b=True,which='major',color='#666666')
			grid(b=True,which='minor',color='#999999')
			savefig(pathname+str(param_array[k])+'.pdf')
			close()

		
		# Plot for all output parameters
		figure()
		flag=0
		for k in range(n_param):
			if param_array[k] not in output_conditions:
				continue
			semilogx(current_array,extracted_matrix[i,:,k],color=colour_dict[param_array[k]],label=param_array[k])
			arrY=output_conditions[param_array[k]]*np.ones(n_current,dtype=float)
			semilogx(current_array,arrY,color=colour_dict[param_array[k]],linestyle='--',label='Output Spec')

			if flag==0:
				max_val=np.max(extracted_matrix[i,:,k])
				min_val=np.min(extracted_matrix[i,:,k])
				flag=1
			else:
				cur_max=np.max(extracted_matrix[i,:,k])
				cur_min=np.min(extracted_matrix[i,:,k])
				if max_val<cur_max:
					max_val=cur_max
				if min_val>cur_min:
					min_val=cur_min

			if output_conditions[param_array[k]]>max_val:
				max_val=output_conditions[param_array[k]]
			if output_conditions[param_array[k]]<min_val:
				min_val=output_conditions[param_array[k]]

		arrX=spec_current*np.ones(100,dtype=float)
		arrY=np.linspace(min_val,max_val,100)
		semilogx(arrX,arrY,color='black',linestyle='--')
		xlabel('Io')
		ylabel(param_array[k])
		legend()
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+'all.pdf')
		close()



	# Second, we will plot Parameters vs Current for all temperatures in the same plot
	multi_colour=['green','blue','yellow','cyan','darkviolet','orange','peru']
	multi_linestyle=['-','--','-.',':']
	
	# Creating the directory
	pathname=file_sub_directory+'all_temp/'
	if not os.path.exists(pathname):
		os.makedirs(pathname)

	for k in range(n_param):
			
		figure()

		for i in range(n_temp):
			
			# Plotting the main parameter
			semilogx(current_array,extracted_matrix[i,:,k],color=multi_colour[i%7],linestyle=multi_linestyle[int(i/7)],label=temp_array[i])
			
			# Calculating the maximum and minimum value of the output
			if i==0:
				max_val=np.max(extracted_matrix[i,:,k])
				min_val=np.min(extracted_matrix[i,:,k])
			else:
				cur_max=np.max(extracted_matrix[i,:,k])
				cur_min=np.min(extracted_matrix[i,:,k])
				if max_val<cur_max:
					max_val=cur_max
				if min_val>cur_min:
					min_val=cur_min
			
		if param_array[k] in output_conditions:
			arrY=output_conditions[param_array[k]]*np.ones(n_current,dtype=float)
			semilogx(current_array,arrY,color='red',label='Output Spec')
			if output_conditions[param_array[k]]>max_val:
				max_val=output_conditions[param_array[k]]
			if output_conditions[param_array[k]]<min_val:
				min_val=output_conditions[param_array[k]]

		# Plotting a vertical line at spec_current
		arrX=spec_current*np.ones(100,dtype=float)
		arrY=np.linspace(min_val,max_val,100)
		semilogx(arrX,arrY,color='black',linestyle='--')

		# Other Plotting Parameters
		xlabel('Io')
		ylabel(param_array[k])
		legend()
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+str(param_array[k])+'.pdf')
		close()
	
	
#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Temperature at different current )
def plot_param_vs_temperature(extracted_matrix,temp_array,current_array,param_array,file_sub_directory,output_conditions,colour_dict,spec_current):
	
	# Finding the length of each array
	n_temp=len(temp_array)
	n_current=len(current_array)
	n_param=len(param_array)

	# First, we will plot Parameters vs Temperature and each plot will contain a single current
	for j in range(n_current):
		
		# Creating the directory
		pathname=file_sub_directory+cf.num_trunc(current_array[j],3)+'/'
		if not os.path.exists(pathname):
			os.makedirs(pathname)

		for k in range(n_param):
			
			# Plot for single parameter
			figure()
			
			# Plotting the main parameter
			plot(temp_array,extracted_matrix[:,j,k],color='g',label=param_array[k])
			
			# Calculating the maximum and minimum value of the output
			max_val=np.max(extracted_matrix[:,j,k])
			min_val=np.min(extracted_matrix[:,j,k])
			if param_array[k] in output_conditions:
				arrY=output_conditions[param_array[k]]*np.ones(n_temp,dtype=float)
				plot(temp_array,arrY,color='r',label='Output Spec')
				if output_conditions[param_array[k]]>max_val:
					max_val=output_conditions[param_array[k]]
				if output_conditions[param_array[k]]<min_val:
					min_val=output_conditions[param_array[k]]

			# Plotting a vertical line at spec_current
			arrX=27*np.ones(100,dtype=float)
			arrY=np.linspace(min_val,max_val,100)
			plot(arrX,arrY,color='black',linestyle='--')

			# Other Plotting Parameters
			xlabel('Temperature')
			ylabel(param_array[k])
			legend()
			grid(b=True,which='major',color='#666666')
			grid(b=True,which='minor',color='#999999')
			savefig(pathname+str(param_array[k])+'.pdf')
			close()

		
		# Plot for all output parameters
		figure()
		flag=0
		for k in range(n_param):
			if param_array[k] not in output_conditions:
				continue
			plot(temp_array,extracted_matrix[:,j,k],color=colour_dict[param_array[k]],label=param_array[k])
			arrY=output_conditions[param_array[k]]*np.ones(n_temp,dtype=float)
			plot(temp_array,arrY,color=colour_dict[param_array[k]],linestyle='--',label='Output Spec')

			if flag==0:
				max_val=np.max(extracted_matrix[:,j,k])
				min_val=np.min(extracted_matrix[:,j,k])
				flag=1
			else:
				cur_max=np.max(extracted_matrix[:,j,k])
				cur_min=np.min(extracted_matrix[:,j,k])
				if max_val<cur_max:
					max_val=cur_max
				if min_val>cur_min:
					min_val=cur_min

			if output_conditions[param_array[k]]>max_val:
				max_val=output_conditions[param_array[k]]
			if output_conditions[param_array[k]]<min_val:
				min_val=output_conditions[param_array[k]]

		arrX=27*np.ones(100,dtype=float)
		arrY=np.linspace(min_val,max_val,100)
		plot(arrX,arrY,color='black',linestyle='--')
		xlabel('Temperature')
		ylabel(param_array[k])
		legend()
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+'all.pdf')
		close()



	# Second, we will plot Parameters vs Current for all temperatures in the same plot
	multi_colour=['green','blue','yellow','cyan','darkviolet','orange','peru']
	multi_linestyle=['-','--','-.',':']
	
	# Creating the directory
	pathname=file_sub_directory+'all_current/'
	if not os.path.exists(pathname):
		os.makedirs(pathname)

	for k in range(n_param):
			
		figure()
		
		for j in range(n_current):
			
			# Plotting the main parameter
			plot(temp_array,extracted_matrix[:,j,k],color=multi_colour[j%7],linestyle=multi_linestyle[int(j/7)],label=cf.num_trunc(current_array[j],3))
			
			# Calculating the maximum and minimum value of the output
			if j==0:
				max_val=np.max(extracted_matrix[:,j,k])
				min_val=np.min(extracted_matrix[:,j,k])
			else:
				cur_max=np.max(extracted_matrix[:,j,k])
				cur_min=np.min(extracted_matrix[:,j,k])
				if max_val<cur_max:
					max_val=cur_max
				if min_val>cur_min:
					min_val=cur_min
			
		if param_array[k] in output_conditions:
			arrY=output_conditions[param_array[k]]*np.ones(n_temp,dtype=float)
			plot(temp_array,arrY,color='red',label='Output Spec')
			if output_conditions[param_array[k]]>max_val:
				max_val=output_conditions[param_array[k]]
			if output_conditions[param_array[k]]<min_val:
				min_val=output_conditions[param_array[k]]

		# Plotting a vertical line at spec_current
		arrX=27*np.ones(100,dtype=float)
		arrY=np.linspace(min_val,max_val,100)
		plot(arrX,arrY,color='black',linestyle='--')

		# Other Plotting Parameters
		xlabel('Temperature')
		ylabel(param_array[k])
		legend()
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+str(param_array[k])+'.pdf')
		close()

#===========================================================================================================================




		
	


