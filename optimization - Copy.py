#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: optimization.py
File Description 	: This file will perform the optimization for different circuit parameters

Functions structure in this file:
	--> save_input_results_optimization
	--> save_output_results_optimization
	--> save_info_single_array_iter
	--> save_info_double_array_iter
	--> save_info
	--> extract_double_array
	--> plot_double_array
	--> extract_single_array
	--> plot_single_array
	--> calc_loss_slope
	--> update_alpha
	--> check_circuit_parameters
	--> update_C2_Rbias
	--> check_stop_alpha
	--> check_stop_loss
	--> moving_avg
	--> opt_single_run
	--> main_opt

	
"""
#===========================================================================================================================
import numpy as np
import common_functions as cf
import optimization_functions_loss as ofl
import optimization_functions_fom as off
import os
from matplotlib import pylab
from pylab import *
#===========================================================================================================================

"""
===========================================================================================================================
-------------------------------------------- Storing Optimization Results -------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that stores input data of the simulation
# Inputs  : optimization_input_parameters
# Outputs : NONE
def save_input_results_optimization(optimization_input_parameters,run_number):

	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/input_data.txt')
	f=open(filename,'a')

	# Storing the results
	if run_number==1:
		f.write('\n\n---------------------- Optimization Parameters -----------------------')
		f.write('\nOptimization Name :'+str(optimization_input_parameters['optimization']['optimization_name']))

		if 'acceptable_solution' in optimization_input_parameters:
			f.write('\n\n---------------------- Acceptable Solution Parameters -----------------------')
			for name in optimization_input_parameters['acceptable_solution']:
				f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['acceptable_solution'][name],3))
	

	f.write('\n\n---------------------- Run Number '+str(run_number)+' -----------------------')
	
	f.write('\nMax Iterations :'+str(optimization_input_parameters['optimization'][run_number]['max_iteration']))
	f.write('\nAlpha Min      :'+str(optimization_input_parameters['optimization'][run_number]['alpha_min']))
	f.write('\nConsec Iter    :'+str(optimization_input_parameters['optimization'][run_number]['consec_iter']))
	f.write('\nAlpha Mult      :'+str(optimization_input_parameters['optimization'][run_number]['alpha_mult']))
	f.write('\nDelta Threshold :'+str(optimization_input_parameters['optimization'][run_number]['delta_threshold']))
	f.write('\nLoss Type       :'+str(optimization_input_parameters['optimization'][run_number]['loss_type']))
	f.write('\nUpdate Check    :'+str(optimization_input_parameters['optimization'][run_number]['update_check']))
	f.write('\nOptimization Type :'+str(optimization_input_parameters['optimization'][run_number]['optimization_type']))

	f.write('\nOptimization Parameters : ')
	for name in optimization_input_parameters['optimization'][run_number]['optimizing_parameters']:
		f.write(str(name)+' ,')

	f.write('\n\n---------------------- Loss Weights -----------------------')
	for name in optimization_input_parameters['optimization'][run_number]['loss_weights']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['optimization'][run_number]['loss_weights'][name],3))

	f.write('\n\n---------------------- Alpha Parameters -----------------------')
	#for name in optimization_input_parameters['optimization'][run_number]['alpha']['values']:
	#	f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['optimization'][run_number]['alpha']['values'][name],3))
	f.write('\nAlpha Value :'+str(optimization_input_parameters['optimization'][run_number]['alpha']['value']))
	f.write('\nAlpha Type  :'+str(optimization_input_parameters['optimization'][run_number]['alpha']['type']))
	f.write('\nAlpha Start :'+str(optimization_input_parameters['optimization'][run_number]['alpha']['start']))
	f.write('\nAlpha End   :'+str(optimization_input_parameters['optimization'][run_number]['alpha']['end']))

	f.write('\n\n---------------------- Simulation Parameters -----------------------')
	for param_name in optimization_input_parameters['optimization']['simulation'][run_number]['standard_parameters']:
		f.write('\n'+param_name+' : '+str(optimization_input_parameters['optimization']['simulation'][run_number]['standard_parameters'][param_name]))
	for param_name in optimization_input_parameters['optimization']['simulation'][run_number]['netlist_parameters']:
		f.write('\n'+param_name+' : '+str(optimization_input_parameters['optimization']['simulation'][run_number]['netlist_parameters'][param_name]))
	
	f.close()


#-----------------------------------------------------------------
# Function that stores output data of the optimization
# Inputs  : optimization_results, optimization_input_parameters
# Outputs : NONE
def save_output_results_optimization(optimization_results,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
	f=open(filename,'a')

	print_dict=optimization_results['optimized_results']
	iter_number=print_dict['iter_number']-1

	run_number=optimization_results['run_number']

	# Storing the results
	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Optimization '+str(run_number)+'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	
	if 'optimization_start' in optimization_results:
		f.write('\n\n--------------------- Optimization Start ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cf.print_output_parameters(f,optimization_results['optimization_start']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cf.print_output_parameters(f,optimization_results['optimization_start']['extracted_parameters'])
	

	f.write('\n-------------------------------------------------------------------')
	if optimization_input_parameters['optimization']['optimization_name']=='loss1':
		f.write('\nMaximum Loss of gain+Io+s11+iip3='+cf.num_trunc(print_dict['loss_max'],3))
		f.write('\nHas the optimization yielded a correct point :'+print_dict['perfect_point'])
		f.write('\nOptimized Point occured at iteration='+str(print_dict['iter_number']))
		f.write('\nOptimized Io Loss='+cf.num_trunc(print_dict['Io_loss'],3))
	
	elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
		f.write('\nMaximum Loss of s11='+cf.num_trunc(print_dict['loss_max'],3))
		f.write('\nOptimized Point occured at iteration='+str(print_dict['iter_number']))
		f.write('\nOptimized FOM in dB='+cf.num_trunc(print_dict['FOM'],3))
	
	f.write('\n\n--------------------- Optimization End ---------------------------------')
	f.write('\n\n---------------- Circuit Parameters ------------------------')
	cf.print_output_parameters(f,optimization_results['circuit_parameters_iter'][iter_number])
	f.write('\n\n---------------- Extracted Parameters ------------------------')
	cf.print_output_parameters(f,optimization_results['extracted_parameters_iter'][iter_number])

	if 'acceptable_solution' in optimization_results:
		f.write('Acceptable Solutions:\n')
		for i in optimization_results['acceptable_solution']:
			f.write(str(i)+' ; ')
	
	f.close()


"""
===========================================================================================================================
-------------------------------------------- Storing Iteration Results ----------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that stores the data of parameters vs iterations in a csv file
# Inputs  : filename_root,filename_name,values_iter,niter
# Outputs : NONE
def save_info_single_array_iter(filename_root,filename_name,values_iter,iter_no):
	
	filename=filename_root+filename_name
	if iter_no==0:
		f=open(filename,'w')
	
	else:
		f=open(filename,'a')

	if iter_no==0:
		f.write('Iteration No,')
		for param in values_iter[0]:
			f.write(str(param)+',')
		f.write('\n')
	
	
	f.write(str(iter_no+1)+',')
	for param in values_iter[iter_no]:
		f.write(str(values_iter[iter_no][param])+',')
	f.write('\n')
	
	f.close()
	
#-----------------------------------------------------------------
# Function that stores the data of loss slopes vs iterations in a csv file
# Inputs  : filename_root,filename_name,values_iter,niter
# Outputs : NONE
def save_info_double_array_iter(filename_root,filename_name,values_iter,iter_no):
	
	if iter_no==0:
		return
	iter_no=iter_no-1
	
	filename=filename_root+filename_name
	
	if iter_no==0:
		f=open(filename,'w')
	
	else:
		f=open(filename,'a')
	
	

	if iter_no==0:
		f.write('Iteration No,')
		for param in values_iter[0]:
			f.write(str(param)+',')
			for categ in values_iter[0][param]:
				f.write(str(categ)+',')
		f.write('\n')
	
	f.write(str(iter_no+1)+',')
	for param in values_iter[iter_no]:
		f.write(str(param)+',')
		for categ in values_iter[iter_no][param]:
			f.write(str(values_iter[iter_no][param][categ])+',')
	f.write('\n')
	
	f.close()

#-----------------------------------------------------------------
# Function that stores the all the simulation data in different csv files
# Inputs  : optimization_input_parameters,optimization_results
# Outputs : NONE
def save_info(optimization_input_parameters,optimization_results,iter_no,flag):
	
	# Creating the folder to store the results
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/Optimization'+str(optimization_results['run_number'])+'/results/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	loss_slope_iter				= optimization_results['loss_slope_iter']
	sensitivity_iter			= optimization_results['sensitivity_iter']
	loss_iter					= optimization_results['loss_iter']
	alpha_parameters_iter		= optimization_results['alpha_parameters_iter']
	extracted_parameters_iter	= optimization_results['extracted_parameters_iter']
	circuit_parameters_iter		= optimization_results['circuit_parameters_iter']
	
	
	save_info_double_array_iter(newpath,'loss_slope.csv',loss_slope_iter,iter_no)
	save_info_double_array_iter(newpath,'sensitivity.csv',sensitivity_iter,iter_no)
	
	if flag==1:
		save_info_single_array_iter(newpath,'loss.csv',loss_iter,iter_no)
		save_info_single_array_iter(newpath,'alpha_parameters.csv',alpha_parameters_iter,iter_no)
		save_info_single_array_iter(newpath,'output_parameters.csv',extracted_parameters_iter,iter_no)
		save_info_single_array_iter(newpath,'circuit_parameters.csv',circuit_parameters_iter,iter_no)
	

"""
===========================================================================================================================
-------------------------------------------- Plotting the Results ---------------------------------------------------------
"""
	
#-----------------------------------------------------------------------------------------------
# Function to extract the data from a 2D dictionary file
def extract_double_array(optimization_results,list_name):
	
	# Assigning the array to store the loss variables
	param_name_list=[]
	variable_name_list=[]

	# Creating variables for the counting no of lines and variables
	n_iter=optimization_results['n_iter']
	
	flag=0
	n_param=0
	for param_name in optimization_results[list_name][1]:
		n_param+=1
		param_name_list.append(param_name)
		if flag==0:
			param_name1=param_name
			flag=1
	
	n_variable=0
	for variable_name in optimization_results[list_name][1][param_name1]:
		n_variable+=1
		variable_name_list.append(variable_name)
	
	
	# Creating array to store the values 
	parameter_matrix=np.zeros((n_iter,n_param,n_variable),dtype=float)
	
	# Extracting the values of the variables
	for i in range(n_iter):
	
		# Storing the variables
		j=0
		for param_name in optimization_results[list_name][i]:
			k=0
			for variable_name in optimization_results[list_name][i][param_name]:
				parameter_matrix[i,j,k]=optimization_results[list_name][i][param_name][variable_name]
		
	return parameter_matrix,param_name_list,variable_name_list

	
#-----------------------------------------------------------------------------------------------
# Plotting 2D array dictionary vs iterations
def plot_double_array(optimization_results,list_name,file_directory):
	
	loss_array,param_name,variable_name=extract_double_array(optimization_results,list_name)
	
	n_iter=loss_array.shape[0]
	n_param=loss_array.shape[1]
	n_variable=loss_array.shape[2]
	
	filename=file_directory+'/plots/'+list_name+'/'
	if not os.path.exists(filename):
		os.makedirs(filename)
	
	# Creating the new arrays
	arrX=np.zeros((n_iter,1),dtype=float)
	
	# Calculating values of new array
	i=0	
	while i<n_iter:
		arrX[i,0]=i+1	
		i+=1
		
	print('Starting Plots '+list_name)
	
	color_dict={0:'r',1:'g',2:'b',3:'c',4:'m',5:'k',6:'y'}
	n_colour=7
	
	for i in range(n_param):
		figure()
		for j in range(n_variable):
			plot(arrX,loss_array[:,i,j],color_dict[int(j%n_colour)],label=variable_name[j])
		xlabel('Iterations')
		ylabel(param_name[i])
		legend()
		grid()
		savefig(filename+param_name[i]+'.pdf')
		close()
	
	print('Plotting Over '+list_name)

	
#-----------------------------------------------------------------------------------------------
# Function to extract the data from 1D dictionary
def extract_single_array(optimization_results,list_name):
	
	# Assigning the array to store the loss variables
	param_name_list=[]
	
	# Creating variables for the counting no of lines and variables
	n_iter=optimization_results['n_iter']
	
	n_param=0
	for param_name in optimization_results[list_name][1]:
		param_name_list.append(param_name)
		n_param+=1
	
	# Creating array to store the values 
	parameter_matrix=np.zeros((n_iter,n_param),dtype=float)
	
	# Extracting the values of the variables
	for i in range(n_iter):
		
		# Storing the variables
		j=0
		for param_name in optimization_results[list_name][i]:
			parameter_matrix[i,j]=optimization_results[list_name][i][param_name]
			j+=1
		
	return parameter_matrix,param_name_list
	
	
#-----------------------------------------------------------------------------------------------
# Plotting parameters vs iterations
def plot_single_array(optimization_results,list_name,file_directory):
	
	loss_array,param_name=extract_single_array(optimization_results,list_name)
	
	n_iter=loss_array.shape[0]
	n_param=loss_array.shape[1]
	
	filename=file_directory+'/plots/'+list_name+'/'
	if not os.path.exists(filename):
    		os.makedirs(filename)
	
	# Creating the new arrays
	arrX=np.zeros((n_iter,1),dtype=float)
	
	# Calculating values of new array
	i=0	
	while i<n_iter:
		arrX[i,0]=i+1	
		i+=1
		
	print('Starting Plots '+list_name)
	
	color_dict={0:'r',1:'g',2:'b',3:'c',4:'m',5:'k',6:'y'}
	n_colour=7
	
	# Figure 1
	figure()
	for i in range(n_param):
		plot(arrX,loss_array[:,i],color_dict[int(i%n_colour)],label=param_name[i])
		annotate(cf.num_trunc(loss_array[n_iter-1,i],3),(arrX[n_iter-1,0],loss_array[n_iter-1,i]))
	xlabel('Iterations')
	ylabel('Parameter')
	legend()
	grid()
	savefig(filename+'all.pdf')
	close()
	
	
	# Figure 2	
	for i in range(n_param):
		figure()
		plot(arrX,loss_array[:,i],'g',label=param_name[i])
		annotate(cf.num_trunc(loss_array[n_iter-1,i],3),(arrX[n_iter-1,0],loss_array[n_iter-1,i]))
		xlabel('Iterations')
		ylabel(param_name[i])
		legend()
		grid()
		savefig(filename+param_name[i]+'.pdf')
		close()
	
	print('Plotting Over '+list_name)


#-----------------------------------------------------------------------------------------------
# Plotting various plots for optimization process
def plot_optimization(optimization_input_parameters,optimization_results,run_number):
	file_directory=optimization_input_parameters['filename']['output']+'/Optimization'+str(run_number)
	plot_single_array(optimization_results,'loss_iter',file_directory)
	plot_single_array(optimization_results,'alpha_parameters_iter',file_directory)
	plot_single_array(optimization_results,'extracted_parameters_iter',file_directory)
	plot_single_array(optimization_results,'circuit_parameters_iter',file_directory)
	
	plot_double_array(optimization_results,'loss_slope_iter',file_directory)
	plot_double_array(optimization_results,'sensitivity_iter',file_directory)


"""
===========================================================================================================================
------------------------------------Defining the functions ----------------------------------------------------------------
"""
	
#-----------------------------------------------------------------------------------------------
# This function updates the values of circuit parameters by trying to minimize loss
# Inputs  : output_conditions,circuit_parameters,loss_dict,extracted_parameters,optimization_input_parameters
# Outputs : circuit_parameters_slope,circuit_parameters_sensitivity
def calc_loss_slope(cir,output_conditions,loss_dict,optimization_input_parameters,run_number):

	loss_weights=optimization_input_parameters['optimization'][run_number]['loss_weights']
	delta_threshold=optimization_input_parameters['optimization'][run_number]['delta_threshold']
	
	# Getting the sensitivity dictionary
	circuit_parameters_sensitivity={}
	for param_name in optimization_input_parameters['optimization'][run_number]['optimizing_parameters']:
		circuit_parameters_sensitivity[param_name]=0
	
	# Creating new dictionaries
	circuit_parameters_initial=cir.circuit_parameters.copy()
	extracted_parameters_initial=cir.extracted_parameters.copy()
	circuit_parameters1=cir.circuit_parameters.copy() # This dictionary will store the values of parameters after increment to calculate the slope
	extracted_parameters1=cir.extracted_parameters.copy() # This dictionary will store the values of parameters after increment to calculate the slope
	circuit_parameters_slope={} # This dictionary will store the values of slope of different losses with change of all circuit parameters
	
	# Calculating the value to update each parameter with
	for param_name in optimization_input_parameters['optimization'][run_number]['optimizing_parameters']:
		
		# Calculating the increment value
		increment_factor=delta_threshold # The value by which parameter increases = increment_factor*parameter
		increment=cir.circuit_parameters[param_name]*increment_factor
	
		# Incrementing the circuit parameter
		circuit_parameters1=cir.circuit_parameters.copy()
		circuit_parameters1[param_name]=circuit_parameters1[param_name]+increment
		
		# Extracting Loss
		cir.update_circuit(circuit_parameters1)
		extracted_parameters1=cir.extracted_parameters.copy()
		
		if optimization_input_parameters['optimization']['optimization_name']=='loss1':
			loss_dict1=ofl.calc_loss_1(extracted_parameters1,output_conditions,loss_weights)
		elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
			loss_dict1=off.calc_fom_1(extracted_parameters1,output_conditions,loss_weights)
		
		# Calculating Slope	
		circuit_parameters_slope[param_name]={}
		for param in loss_dict:
			circuit_parameters_slope[param_name][param]=(loss_dict1[param]-loss_dict[param])/increment
			
		circuit_parameters_sensitivity[param_name]={}
		
		# Calculating Sensitivity
		for categ_name in optimization_input_parameters['optimization'][run_number]['output_parameters_list']:
			initial_param=extracted_parameters_initial[categ_name]
			final_param=extracted_parameters1[categ_name]
			percent_change=(final_param-initial_param)/(initial_param*increment_factor)
			circuit_parameters_sensitivity[param_name][categ_name]=percent_change
	
	cir.circuit_parameters=circuit_parameters_initial.copy()
	cir.extracted_parameters=extracted_parameters_initial.copy()
		
	return circuit_parameters_slope,circuit_parameters_sensitivity
	
#-----------------------------------------------------------------------------------------------
# This function updates the value of alpha after each iteration
# Inputs  : loss_iter,alpha,i,alpha_mult,optimization_type,optimization_input_parameters
# Outputs : alpha
def update_alpha(loss_iter,alpha,i,alpha_mult,optimization_type,optimization_input_parameters,run_number):

	n_iter=optimization_input_parameters['optimization'][run_number]['max_iteration']-1
	alpha_start=optimization_input_parameters['optimization'][run_number]['alpha']['start']
	alpha_end=optimization_input_parameters['optimization'][run_number]['alpha']['end']

	if optimization_input_parameters['optimization'][run_number]['alpha']['type']=='Linear':
		alpha=alpha_start+((alpha_end-alpha_start)*(i+1)/n_iter)
		print(alpha)

	elif optimization_input_parameters['optimization'][run_number]['alpha']['type']=='Log':
		alpha_start_log=np.log(alpha_start)
		alpha_end_log=np.log(alpha_end)
		alpha_log=alpha_start_log+((alpha_end_log-alpha_start_log)*(i+1)/n_iter)
		alpha=np.exp(alpha_log)
		
	else:
		# Checking criteria for reducing threshold
		if i>0:
			if loss_iter[i]['loss']>=loss_iter[i-1]['loss'] and optimization_type==0:
				alpha*=alpha_mult
			elif loss_iter[i]['loss']<=loss_iter[i-1]['loss'] and optimization_type==1:
				alpha*=alpha_mult
		
	return alpha

#-----------------------------------------------------------------------------------------------
# This function updates circuit parameters to previous circuit parameters if loss increases
# Inputs  : old_circuit_parameters,circuit_parameters,loss_iter,update_check,i,optimization_type
# Outputs : circuit_parameters, old_circuit_parameters
def check_circuit_parameters(old_circuit_parameters,cir,loss_iter,update_check,i,optimization_type):

	# Checking criteria for reducing threshold
	if update_check==1:
		if i>0:
			if loss_iter[i]['loss']>=loss_iter[i-1]['loss'] and optimization_type==0:
				cir.circuit_parameters=old_circuit_parameters.copy()
			elif loss_iter[i]['loss']<=loss_iter[i-1]['loss'] and optimization_type==1:
				cir.circuit_parameters=old_circuit_parameters.copy()
	return cir.circuit_parameters

#-----------------------------------------------------------------------------------------------
# Checking stopping condition ( if alpha<alpha_min )
# Inputs  : loss_iter,alpha,i,alpha_min
# Outputs : 1 if we need to stop iterations
def check_stop_alpha(loss_iter,alpha,i,alpha_min):

	if alpha_min<0:
		return 0
	if i>1:
		if alpha<=alpha_min:
			return 1
	return 0

#-----------------------------------------------------------------------------------------------
# Checking stopping condition ( if loss increases for n_iter number of iterations )
# Inputs  : loss_iter,i,n_iter,optimization_type
# Outputs : 1 if we need to stop iterations
def check_stop_loss(loss_iter,i,n_iter,optimization_type):

	flag=0
	if n_iter<1:
		return 0
	if i>n_iter:
		flag=1
		for j in range(n_iter):
			if loss_iter[i-j]['loss']<loss_iter[i-j-1]['loss'] and optimization_type==0:
				flag=0
			elif loss_iter[i-j]['loss']>loss_iter[i-j-1]['loss'] and optimization_type==1:
				flag=0
	return flag
	

"""
===========================================================================================================================
------------------------------------------- Output Functions --------------------------------------------------------------
"""
		
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to do optimization for a single run 
# Inputs  : circuit_parameters,extracted_parameters,optimization_input_parameters,run_number
# Outputs : circuit_parameters,extracted_parameters
def opt_single_run(cir,optimization_input_parameters,run_number):

	optimization_results={}
	optimization_results['run_number']=run_number

	# Defining some values
	i=0
	
	output_conditions=optimization_input_parameters['output_conditions']
	
	loss_weights      = optimization_input_parameters['optimization'][run_number]['loss_weights']
	alpha_min         = optimization_input_parameters['optimization'][run_number]['alpha_min']
	consec_iter       = optimization_input_parameters['optimization'][run_number]['consec_iter']
	alpha_mult        = optimization_input_parameters['optimization'][run_number]['alpha_mult']
	max_iteration     = optimization_input_parameters['optimization'][run_number]['max_iteration']
	delta_threshold   = optimization_input_parameters['optimization'][run_number]['delta_threshold']
	loss_type         = optimization_input_parameters['optimization'][run_number]['loss_type']
	optimization_type = optimization_input_parameters['optimization'][run_number]['optimization_type']
	
	#alpha_parameters         = optimization_input_parameters['optimization'][run_number]['alpha']['values']
	#alpha_parameters_initial = optimization_input_parameters['optimization'][run_number]['alpha']['values'].copy()

	alpha         = optimization_input_parameters['optimization'][run_number]['alpha']['value']
	alpha_initial = optimization_input_parameters['optimization'][run_number]['alpha']['value']
	
	# Creating old circuit parameters
	old_circuit_parameters=cir.circuit_parameters.copy() # This dictionary will store the value of parameters for previous iterations
	
	# Creating the dictionaries
	loss_iter={} 				# This dictionary will store the value of all loss values for different iterations
	loss_slope_iter={} 			# This dictionary will store the value of slope of losses for all parameters for different iterations
	alpha_parameters_iter={} 	# This dictionary will store the value of threshold for different iterations
	extracted_parameters_iter={}# This dictionary will store the value of output parameters for different iterations
	circuit_parameters_iter={} 	# This dictionary will store the value of circuit parameters for different iterations
	sensitivity_iter={}			# This dictionary will store the value of output parameter sensitivty for different circuit parameters 
	check_loss=1	
	
	# Running Eldo
	cir.run_circuit()
	
	# Storing the Circuit and Extracted Parameters
	optimization_results['optimization_start']={}
	optimization_results['optimization_start']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['optimization_start']['extracted_parameters']=cir.extracted_parameters.copy()

	# Calculating loss
	if optimization_input_parameters['optimization']['optimization_name']=='loss1':
		loss_iter[-1]=ofl.calc_loss_1(cir.extracted_parameters,output_conditions,loss_weights)
	elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
		loss_iter[-1]=off.calc_fom_1(cir.extracted_parameters,output_conditions,loss_weights)
	
	# Printing the values of loss before optimization
	print('-----------------------------Before Iteration---------------------------------')
	cf.print_loss_parameters(loss_iter[-1])


	# Storing optimization results
	optimization_results['loss_iter']					= loss_iter
	optimization_results['loss_slope_iter']				= loss_slope_iter
	optimization_results['alpha_parameters_iter']		= alpha_parameters_iter
	optimization_results['extracted_parameters_iter']	= extracted_parameters_iter
	optimization_results['circuit_parameters_iter']		= circuit_parameters_iter
	optimization_results['sensitivity_iter']			= sensitivity_iter


	#--------------------------- Performing the iterations ---------------------------------
	while i<max_iteration:
	
		# Checking if there is extra loss from output conditions
		if optimization_input_parameters['optimization']['optimization_name']=='loss1':
			check_loss=ofl.calc_check_loss(loss_iter,i,loss_type)
		elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
			check_loss=off.calc_check_loss(loss_iter,i,loss_type)
		
		
		# Calculating the slope of loss and output sensitivity and updating the circuit parameters
		circuit_parameters_slope,circuit_parameters_sensitivity=calc_loss_slope(cir,output_conditions,loss_iter[i-1],optimization_input_parameters,run_number)
		if optimization_input_parameters['optimization']['optimization_name']=='loss1':
			ofl.update_circuit_parameters(cir,circuit_parameters_slope,check_loss,optimization_input_parameters,run_number)
		elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
			off.update_circuit_parameters(cir,circuit_parameters_slope,check_loss,optimization_input_parameters,run_number)
		

		# Extracting output parameters for new circuit parameters
		cir.run_circuit()


		# Updating different dictionaries
		if optimization_input_parameters['optimization']['optimization_name']=='loss1':
			loss_iter[i]=ofl.calc_loss_1(cir.extracted_parameters,output_conditions,loss_weights)
		elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
			loss_iter[i]=off.calc_fom_1(cir.extracted_parameters,output_conditions,loss_weights)
			

		# Storing some parameters
		#alpha_parameters_iter[i]	= alpha_parameters.copy()
		alpha_parameters_iter[i]	= {'alpha':alpha}
		loss_slope_iter[i-1]		= circuit_parameters_slope.copy()
		sensitivity_iter[i-1]		= circuit_parameters_sensitivity.copy()
		circuit_parameters_iter[i]	= cir.circuit_parameters.copy()
		extracted_parameters_iter[i]= cir.extracted_parameters.copy()

		# Saving Results of Each Iteration
		save_info(optimization_input_parameters,optimization_results,i,1)
		

		# Opening the Run_Status File
		f=open(optimization_input_parameters['filename']['run_status'],'a')
		f.write('Iteration Number:'+str(i+1)+'\n')
		f.close()


		# Printing the values of loss for given iteration
		print('\n-----------------------------Iteration Number ',i+1,'-------------------------------')
		cf.print_loss_parameters(loss_iter[i])
		

		# Updating the value of alpha	
		#alpha_parameters['common']=update_alpha(loss_iter,alpha_parameters['common'],i,alpha_mult,optimization_type,optimization_input_parameters,run_number)
		alpha=update_alpha(loss_iter,alpha,i,alpha_mult,optimization_type,optimization_input_parameters,run_number)
		

		# Updating the value of circuit_parameters based on loss increase
		old_circuit_parameters=check_circuit_parameters(old_circuit_parameters,cir,loss_iter,optimization_input_parameters['optimization'][run_number]['update_check'],i,optimization_type)


		# Checking for stopping condition
		flag_alpha=check_stop_alpha(loss_iter,alpha,i,alpha_min)
		flag_loss=check_stop_loss(loss_iter,i,consec_iter,optimization_type)
		

		# Incrementing i and breaking the loop if necessary
		i+=1
		if flag_loss==1 or flag_alpha==1:
			break
	

	# Calculating slope and sensitivity
	circuit_parameters_slope,circuit_parameters_sensitivity=calc_loss_slope(cir,output_conditions,loss_iter[i-1],optimization_input_parameters,run_number)
	loss_slope_iter[i-1]=circuit_parameters_slope.copy()
	sensitivity_iter[i-1]=circuit_parameters_sensitivity.copy()

	
	# Storing the final results
	optimization_results['n_iter']=i
	save_info(optimization_input_parameters,optimization_results,i,0)
	

	# Resetting the value of alpha
	#optimization_input_parameters['optimization'][run_number]['alpha']['values']=alpha_parameters_initial.copy()
	optimization_input_parameters['optimization'][run_number]['alpha']['value']=alpha_initial


	# Finding the best optimization results
	if optimization_input_parameters['optimization']['optimization_name']=='loss1':
		optimization_results['optimized_results']=ofl.check_best_solution(optimization_results,0)
	elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
		optimization_results['optimized_results']=off.check_best_solution(optimization_results,0)
		optimization_results['acceptable_solution']=off.check_acceptable_solutions(optimization_results,optimization_input_parameters)


	# Assigning circuit_parameters and extracted_parameters with the best result
	print_dict=optimization_results['optimized_results']
	iter_number=print_dict['iter_number']-1

	cir.update_circuit(optimization_results['circuit_parameters_iter'][iter_number].copy())
	
	# Printing the values
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)

	# Storing the results of Optimization
	save_output_results_optimization(optimization_results,optimization_input_parameters)

	# Plotting the results of optimization
	plot_optimization(optimization_input_parameters,optimization_results,run_number)
	

	

#---------------------------------------------------------------------------------------------------------------------------
# Function to do optimization for multiple runs
# Inputs  : circuit_parameters,extracted_parameters,optimization_input_parameters,timing_results
# Outputs : circuit_parameters,extracted_parameters
def main_opt(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['optimization']['run']=='NO':
		return

	n_runs=optimization_input_parameters['optimization']['n_runs']
	
	# Storing the starting time
	timing_results['optimization']={}
	timing_results['optimization']['overall']={}
	timing_results['optimization']['overall']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Main Optimization ******************************************************')

	for i in range(1,1+n_runs):

		save_input_results_optimization(optimization_input_parameters,i)
		
		# Opening the Run_Status File
		f=open(optimization_input_parameters['filename']['run_status'],'a')
		f.write('Optimization '+str(i)+' Start\n Time : '+str(datetime.datetime.now())+'\n\n')
		f.close()

		# Storing the starting time
		timing_results['optimization'][i]={}
		timing_results['optimization'][i]['start']=datetime.datetime.now()

		cir.update_simulation_parameters(optimization_input_parameters['optimization']['simulation'][i])
		opt_single_run(cir,optimization_input_parameters,i)

		# Storing the optimization completion time
		timing_results['optimization'][i]['stop']=datetime.datetime.now()

		# Opening the Run_Status File
		f=open(optimization_input_parameters['filename']['run_status'],'a')
		f.write('Optimization '+str(i)+' End\n Time : '+str(datetime.datetime.now())+'\n\n')
		f.close()

	timing_results['optimization']['overall']['stop']=datetime.datetime.now()

#===========================================================================================================================
