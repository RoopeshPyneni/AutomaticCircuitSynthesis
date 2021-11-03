#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: optimization_functions_loss.py
File Description 	: This file will contain the functions for Loss Optimization

Functions structure in this file:
	--> ramp_func
	--> calc_fom_1
	--> update_circuit_parameters
	--> calc_check_loss
	--> check_best_solution
	
"""
#===========================================================================================================================


#===========================================================================================================================
#------------------------------------Defining the functions -----------------------------------------

#-----------------------------------------------------------------------------------------------
# This is the ramp function
# Inputs  : x
# Outputs : r(x)
def ramp_func(x):
	if x>0:
		return x
	else:
		return 0
	
#===========================================================================================================================
#--------------------------------------------Output Functions---------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# This function calculates the loss for Io Optimization
# Inputs  : extracted_parameters,output_conditions,loss_weights
# Outputs : loss_dict
def calc_loss_1(extracted_parameters,output_conditions,loss_weights):
	
	# Extracted Values
	gain=extracted_parameters['gain_db']
	#iip3=extracted_parameters['iip3_dbm']
	s11=extracted_parameters['s11_db']
	nf=extracted_parameters['nf_db']
	Io=extracted_parameters['Io']
	
	# Reference Values
	gain_ref=output_conditions['gain_db']
	#iip3_ref=output_conditions['iip3_dbm']
	s11_ref=output_conditions['s11_db']
	nf_ref=output_conditions['nf_db']
	
	#Defining the weights to calculate Loss
	A1=loss_weights['gain_db']	# Weight for gain
	#A2=loss_weights['iip3_dbm']	# Weight for iip3
	A3=loss_weights['s11_db']	# Weight for s11
	A4=loss_weights['nf_db']	# Weight for nf
	A5=loss_weights['Io']	# Weight for Io
	
	# Calculating Loss
	loss_gain=A1*ramp_func(gain_ref-gain)
	#loss_iip3=A2*ramp_func(iip3_ref-iip3)
	loss_s11=A3*ramp_func(s11-s11_ref)
	loss_nf=A4*ramp_func(nf-nf_ref)
	loss_Io=A5*Io
	#loss=loss_gain+loss_iip3+loss_s11+loss_nf+loss_Io
	loss=loss_gain+loss_s11+loss_nf+loss_Io
	#loss_dict={'loss':loss,'loss_gain':loss_gain,'loss_iip3':loss_iip3,'loss_s11':loss_s11,'loss_nf':loss_nf,'loss_Io':loss_Io}
	loss_dict={'loss':loss,'loss_gain':loss_gain,'loss_s11':loss_s11,'loss_nf':loss_nf,'loss_Io':loss_Io}
	
	return loss_dict

#-----------------------------------------------------------------------------------------------
# This function updates the values of circuit parameters by trying to minimize loss
# Inputs  : circuit_parameters,circuit_parameters_slope,check_loss,optimization_input_parameters
# Outputs : circuit_parameters
def update_circuit_parameters(cir,circuit_parameters_slope,check_loss,optimization_input_parameters):

	alpha_parameters=optimization_input_parameters['optimization']['alpha']['values']

	# Calculating the value to update each parameter with
	for param_name in circuit_parameters_slope:
		
		# Calculating the Increment Value
		if check_loss==-1:
			change=circuit_parameters_slope[param_name]['loss']*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
		elif check_loss==1:
			change=circuit_parameters_slope[param_name]['loss_Io']*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
		else:
			change=(circuit_parameters_slope[param_name]['loss']-circuit_parameters_slope[param_name]['loss_Io'])
			change=change*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
	
	
		# Checking if the parameter is updated by a large value
		change_limit=0.25 # If the incremented value is more than +- change_limit*parameter_name, then we will limit the change
		if change>change_limit*cir.circuit_parameters[param_name]:
			change=change_limit*cir.circuit_parameters[param_name]
		if change<-1*change_limit*cir.circuit_parameters[param_name]:
			change=-1*change_limit*cir.circuit_parameters[param_name]
		
		
		# Updating circuit_parameters
		cir.circuit_parameters[param_name]=cir.circuit_parameters[param_name]-change
		
	
#-----------------------------------------------------------------------------------------------
# This function will check the loss of gain, iip3, nf, and s11
# Inputs  : loss_iter,i,loss_type
# Outputs : check_loss ( -1 if s11 loss is 0 )
def calc_check_loss(loss_iter,i,loss_type):

	if loss_type==0:
		if loss_iter[i-1]['loss']==loss_iter[i-1]['loss_Io']:
			check_loss=1
		else:
			check_loss=0
				
	elif loss_type==1:
		check_loss=-1
		
	return check_loss
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to check the best solution
# Inputs  : optimization_results,loss_max
# Outputs : opt_dict
def check_best_solution(optimization_results,loss_max):

	# Defining some values
	n_iter=optimization_results['n_iter']
	iter_min=0
	loss_Io_min=optimization_results['loss_iter'][0]['loss_Io']

	if (optimization_results['loss_iter'][0]['loss']-optimization_results['loss_iter'][0]['loss_Io'])>loss_max:
		flag=0
	else:
		flag=1
	
	for i in range(1,n_iter):
		if (optimization_results['loss_iter'][i]['loss']-optimization_results['loss_iter'][i]['loss_Io'])>loss_max:
			continue

		if flag==0 or (flag==1 and optimization_results['loss_iter'][i]['loss_Io']<loss_Io_min):
			iter_min=i
			loss_Io_min=optimization_results['loss_iter'][i]['loss_Io']
			flag=1

	# Creating output dictionary
	opt_dict={}
	opt_dict['loss_max']=loss_max
	opt_dict['iter_number']=iter_min+1
	opt_dict['Io_loss']=loss_Io_min
	
	return opt_dict
	

#===========================================================================================================================


		
	


