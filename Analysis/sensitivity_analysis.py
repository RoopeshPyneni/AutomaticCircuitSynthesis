#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the sensitivity analysis by varying each circuit parameter by 1%
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
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
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
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(param,extracted_parameters_sensitivity,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/extracted_parameters_sensitivity.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(param))
	for param_name in extracted_parameters_sensitivity:
		f.write(','+str(extracted_parameters_sensitivity[param_name]))	# Writing the values to the file
	f.write('\n')
	f.close()
	


"""
===========================================================================================================================
--------------------------------------------Output Functions---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the temperature analysis
# Input: circuit_parameters, extracted_parameters, optimization_input_parameters
# Output: circuit_parameters, extracted_parameters
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
	
	initial_extracted_parameters=cir.extracted_parameters.copy()
	initial_circuit_parameters=cir.circuit_parameters.copy()

	# Creating Dictionaries to Store Values
	extracted_parameters_iter={}
	
	# Writing the values to output files
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(cir.extracted_parameters,optimization_input_parameters)
	
	# Performing the analysis
	for param in initial_circuit_parameters:
		current_circuit_parameters=initial_circuit_parameters.copy()
		current_circuit_parameters[param]*=1.01
		cir.update_circuit(current_circuit_parameters)
		extracted_parameters_sensitivity={}
		for ext_param in cir.extracted_parameters:
			if initial_extracted_parameters[ext_param]==0:
				extracted_parameters_sensitivity[ext_param]=0
			else:
				extracted_parameters_sensitivity[ext_param]=100*(cir.extracted_parameters[ext_param]-initial_extracted_parameters[ext_param])/(initial_extracted_parameters[ext_param])
		update_extracted_parameters(param,extracted_parameters_sensitivity,optimization_input_parameters)		# Writing the values to the output file
		
	# Restoring the value of initial extracted and circuit parameters
	cir.update_circuit(initial_circuit_parameters.copy())

	# Plotting the graphs
	#file_directory=optimization_input_parameters['filename']['output']
	#spec_current=cir.circuit_parameters['Io']
	#plot_sensitivity_analysis(extracted_parameters_iter,file_directory,spec_current)

	# Storing the starting time
	timing_results['sensitivity_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Sensitivity Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	
	
#===========================================================================================================================


"""
====================================================================================================================================================================
"""

"""
#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
# Inputs  : extracted_parameters_iter
# Outputs : extracted_matrix, temp_array, current_array, param_array
def plot_sensitivity_analysis(extracted_parameters_iter,file_directory,spec_current):
	
	file_sub_directory=file_directory+'/Sensitivity_Analysis/Plots/'
	extracted_matrix,temp_array,current_array,param_array=extract_temp_analysis(extracted_parameters_iter)

	output_conditions={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}
	colour_dict={'iip3_dbm':'r','nf_db':'g','s11_db':'b','gain_db':'m'}

	plot_param_vs_current(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_current/',output_conditions,colour_dict,spec_current)
	plot_param_vs_temperature(extracted_matrix,temp_array,current_array,param_array,file_sub_directory+'X_temperature/',output_conditions,colour_dict,spec_current)

#-----------------------------------------------------------------------------------------------
# Function to extract the data from extracted_parameters_iter dictionary and store it in the form of a matrix
# Inputs  : extracted_parameters_iter
# Outputs : extracted_matrix, temp_array, current_array, param_array
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
# Inputs  : extracted_matrix, temp_array, current_array, param_array, file_sub_directory, output_conditions, colour_dict, spec_current
# Outputs : NONE
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
# Inputs  : extracted_matrix, temp_array, current_array, param_array, file_sub_directory, output_conditions, colour_dict, spec_current
# Outputs : NONE
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
"""
#===========================================================================================================================
