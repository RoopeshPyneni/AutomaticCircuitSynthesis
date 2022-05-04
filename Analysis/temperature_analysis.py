#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the temperature analysis by sweeping Temperature and Io
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
	
	newpath =filename+'/Temperature_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Temperature_Analysis/Results/initial_circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in initial_circuit_parameters:
		f.write(str(param_name)+'\t'+str(initial_circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
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
def write_extracted_parameters_initial(extracted_parameters,optimization_input_parameters,temp):
	
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis/Results/extracted_parameters_'+str(temp)+'.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Current Mirror Io,')
	for param_name in extracted_parameters:
		f.write(param_name+',')	# Writing the names in the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
def update_extracted_parameters(extracted_parameters,optimization_input_parameters,temp,current):
	filename=optimization_input_parameters['filename']['output']+'/Temperature_Analysis/Results/extracted_parameters_'+str(temp)+'.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(current)+',')
	for param_name in extracted_parameters:
		f.write(str(extracted_parameters[param_name])+',')	# Writing the values to the file
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Get the required data from the extracted parameters
def get_single_temp_extracted_parameters(extracted_parameters,temp):
	
	final_extracted_parameters={}
	for param in extracted_parameters:
		str_temp=str(temp)
		len_temp=len(str_temp)
		if str_temp==param[:len_temp]:
			final_extracted_parameters[param[len_temp+1:]]=extracted_parameters[param]
	
	return final_extracted_parameters
	

"""
============================================================================================================================
------------------------------------- Output Functions ---------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the temperature analysis
def temperature_analysis(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['temperature_analysis']['run']=='NO':
		return

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['temperature_analysis']={}
	timing_results['temperature_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Temperature Analysis ***************************************************')
	
	cir.update_simulation_parameters(optimization_input_parameters['temperature_analysis']['simulation'])
	
	initial_circuit_parameters_initial=cir.get_initial_circuit_parameters()
	circuit_parameters_initial=cir.get_circuit_parameters()
	extracted_parameters_initial=cir.get_extracted_parameters()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Creating an array for Io sweep
	room_temp_current_log=np.log10(cir.initial_circuit_parameters['Io'])
	start_current=np.log10(optimization_input_parameters['temperature_analysis']['start_current'])
	stop_current=np.log10(optimization_input_parameters['temperature_analysis']['stop_current'])
	n_current=optimization_input_parameters['temperature_analysis']['n_current']
	if n_current==1:
		current_array=np.array([cir.initial_circuit_parameters['Io']])
	else:
		current_array=np.logspace(room_temp_current_log+start_current,room_temp_current_log+stop_current,n_current)
	
	print('\n\n----- CURRENT ARRAY -------')
	
	# Writing the values to output files
	write_initial_circuit_parameters(cir.initial_circuit_parameters,optimization_input_parameters)
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	
	# Performing the analysis
	temp_array=cir.circuit_initialization_parameters['simulation']['standard_parameters']['temp_list']

	flag=0
	for current in current_array:
		cir.initial_circuit_parameters['Io']=current
		cir.run_circuit()

		if flag==0:
			flag=1
			for temp in temp_array:
				write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters,temp)
		
		for temp in temp_array:
			final_extracted_parameters=get_single_temp_extracted_parameters(cir.extracted_parameters,temp)
			update_extracted_parameters(final_extracted_parameters,optimization_input_parameters,temp,current)		# Writing the values to the output file
			if temp not in extracted_parameters_iter:
				extracted_parameters_iter[temp]={}
			extracted_parameters_iter[temp][current]=final_extracted_parameters.copy()

	# Restoring the value of parameters
	cir.update_circuit_state(initial_circuit_parameters_initial,circuit_parameters_initial,extracted_parameters_initial)

	# Plotting the graphs
	file_directory=optimization_input_parameters['filename']['output']
	spec_current=cir.initial_circuit_parameters['Io']
	plot_temp_analysis(extracted_parameters_iter,file_directory,spec_current)

	# Storing the starting time
	timing_results['temperature_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Temperature Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	

"""
============================================================================================================================
------------------------------------- Plotting Functions -------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
def plot_temp_analysis(extracted_parameters_iter,file_directory,spec_current):
	
	file_sub_directory=file_directory+'/Temperature_Analysis/Plots/'
	extracted_matrix,temp_array,current_array,param_array=extract_temp_analysis(extracted_parameters_iter)

	output_conditions={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}
	colour_dict={'iip3_dbm':'r','nf_db':'g','s11_db':'b','gain_db':'m'}

	plot_param_vs_current(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_current/',output_conditions,colour_dict,spec_current)
	plot_param_vs_temperature(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_temperature/',output_conditions,colour_dict,spec_current)

#-----------------------------------------------------------------------------------------------
# Function to extract the data from extracted_parameters_iter dictionary and store it in the form of a matrix
def extract_temp_analysis(extracted_parameters_iter):
	
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
# Plotting results ( Parameters vs current at different temperatures )
def plot_param_vs_current(extracted_matrix,temp_array,current_array,param_array,file_sub_directory,output_conditions,colour_dict,spec_current):
	
	# Finding the length of each array
	n_temp=len(temp_array)
	n_current=len(current_array)
	n_param=len(param_array)

	# Finding the limits of X-axis
	percent_cover=0.7
	extra_addition=1.1
	X_arr_start=current_array[0]
	X_arr_stop=current_array[n_current-1]
	X_start=X_arr_start/extra_addition
	log_diff=np.log10(X_arr_stop)-np.log10(X_arr_start)
	X_stop=extra_addition*(10**(np.log10(X_arr_start)+(log_diff/percent_cover)))

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
			xlim([X_start,X_stop])
			legend(loc='upper right', bbox_to_anchor=(1,1.1))
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
		xlim([X_start,X_stop])
		legend(loc='upper right', bbox_to_anchor=(1,1.1))
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+'all.pdf')
		close()



	# Second, we will plot Parameters vs Current for all temperatures in the same plot
	multi_colour=['green','blue','lime','cyan','darkviolet','orange','peru']
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
		#xlim([X_start,X_stop])
		legend(loc='upper right', bbox_to_anchor=(1,1.1))
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

	# Finding the limits of X-axis
	percent_cover=0.7
	extra_addition=10
	X_arr_start=temp_array[0]
	X_arr_stop=temp_array[n_temp-1]
	X_start=X_arr_start-extra_addition
	X_stop=extra_addition+X_arr_start+((X_arr_stop-X_arr_start)/percent_cover)

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
			#xlim([X_start,X_stop])
			legend(loc='upper right', bbox_to_anchor=(1,1.1))
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
		xlim([X_start,X_stop])
		legend(loc='upper right', bbox_to_anchor=(1,1.1))
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+'all.pdf')
		close()



	# Second, we will plot Parameters vs Current for all temperatures in the same plot
	multi_colour=['green','blue','lime','cyan','darkviolet','orange','peru']
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
		xlim([X_start,X_stop])
		legend(loc='upper right', bbox_to_anchor=(1,1.1))
		grid(b=True,which='major',color='#666666')
		grid(b=True,which='minor',color='#999999')
		savefig(pathname+str(param_array[k])+'.pdf')
		close()


#===========================================================================================================================
