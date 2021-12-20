#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will contain the functions for FOM Optimization
"""

#===========================================================================================================================
import common_functions as cf

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
# This function calculates the loss for FoM
# Inputs  : extracted_parameters,output_conditions,loss_weights
# Outputs : fom_dict
def calc_fom_1(extracted_parameters,output_conditions,loss_weights):
	
	# Extracted Values
	gain_db=extracted_parameters['gain_db']
	iip3_dbm=extracted_parameters['iip3_dbm']
	s11_db=extracted_parameters['s11_db']
	nf_db=extracted_parameters['nf_db']
	P=extracted_parameters['p_source']
	freq=extracted_parameters['freq']/1e9
	
	nf=cf.db_to_normal(nf_db)
	
	# Calculating FoM
	fom_gain=gain_db
	fom_iip3=iip3_dbm
	fom_freq=cf.normal_to_db(freq)
	fom_nf=-1*cf.normal_to_db(nf-1)
	fom_P=-1*cf.normal_to_dbm(P)
	fom_s11=-1*loss_weights['s11_db']*ramp_func(s11_db-output_conditions['s11_db'])
	fom=fom_gain+fom_iip3+fom_freq+fom_nf+fom_P+fom_s11
	
	fom_dict={'loss':fom,'loss_gain':fom_gain,'loss_freq':fom_freq,'loss_nf':fom_nf,'loss_iip3':fom_iip3,'loss_P':fom_P,'loss_s11':fom_s11}
	
	return fom_dict

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
			change=circuit_parameters_slope[param_name]['loss']
			change=change*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
		elif check_loss==1:
			change=circuit_parameters_slope[param_name]['loss']-circuit_parameters_slope[param_name]['loss_s11']
			change=change*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
		else:
			change=(circuit_parameters_slope[param_name]['loss_s11'])
			change=change*(cir.circuit_parameters[param_name]**2)*alpha_parameters['common']*alpha_parameters[param_name]
	
	
		# Checking if the parameter is updated by a large value
		change_limit=0.25 # If the incremented value is more than +- change_limit*parameter_name, then we will limit the change
		if change>change_limit*cir.circuit_parameters[param_name]:
			change=change_limit*cir.circuit_parameters[param_name]
		if change<-1*change_limit*cir.circuit_parameters[param_name]:
			change=-1*change_limit*cir.circuit_parameters[param_name]
		
		
		# Updating circuit_parameters
		cir.circuit_parameters[param_name]=cir.circuit_parameters[param_name]+change
		

#-----------------------------------------------------------------------------------------------
# This function will check the loss of s_11
# Inputs  : loss_iter,i,loss_type
# Outputs : check_loss ( -1 if s11 loss is 0 )
def calc_check_loss(loss_iter,i,loss_type):

	if loss_type==0:
		if loss_iter[i-1]['loss_s11']==0:
			check_loss=-1
		else:
			check_loss=0
				
	elif loss_type==1:
		check_loss=-1
		
	return check_loss	
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to check the best solution
# Inputs  : optimization_results,loss_s11_max
# Outputs : opt_dict
def check_best_solution(optimization_results,loss_s11_max):

	# Defining some values
	n_iter=optimization_results['n_iter']
	iter_min=0
	loss_max=optimization_results['loss_iter'][0]['loss']-optimization_results['loss_iter'][0]['loss_s11']

	if optimization_results['loss_iter'][0]['loss_s11']<loss_s11_max:
		flag=0
	else:
		flag=1
	
	for i in range(1,n_iter):
		if optimization_results['loss_iter'][i]['loss_s11']<loss_s11_max:
			continue

		if flag==0 or (flag==1 and (optimization_results['loss_iter'][i]['loss']-optimization_results['loss_iter'][i]['loss_s11'])>loss_max):
			iter_min=i
			loss_max=(optimization_results['loss_iter'][i]['loss']-optimization_results['loss_iter'][i]['loss_s11'])
			flag=1

	# Creating output dictionary
	opt_dict={}
	opt_dict['loss_max']=loss_s11_max
	opt_dict['iter_number']=iter_min+1
	opt_dict['FOM']=loss_max
	
	return opt_dict
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to check the acceptable solutions ( those that satisfy a given criteria )
# Inputs  : optimization_results,optimization_input_parameters
# Outputs : acceptable_iter
def check_acceptable_solutions(optimization_results,optimization_input_parameters):

	# Defining some values
	n_iter=optimization_results['n_iter']
	acceptable_iter=[]

	for i in range(n_iter):
		if optimization_results['output_parameters_iter'][i]['s11_db']>optimization_input_parameters['acceptable_solution']['s11_db']:
			continue	
		if optimization_results['output_parameters_iter'][i]['gain_db']<optimization_input_parameters['acceptable_solution']['gain_db']:
			continue
		if optimization_results['output_parameters_iter'][i]['iip3_dbm']<optimization_input_parameters['acceptable_solution']['iip3_dbm']:
			continue
		if optimization_results['output_parameters_iter'][i]['nf_db']>optimization_input_parameters['acceptable_solution']['nf_db']:
			continue
		if optimization_results['output_parameters_iter'][i]['p_source']>optimization_input_parameters['acceptable_solution']['p_source']:
			continue
		acceptable_iter.append(i)
		
	return acceptable_iter
	

#===========================================================================================================================


		
	


