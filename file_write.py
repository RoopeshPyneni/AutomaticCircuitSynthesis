#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Writing File:
"""
#===========================================================================================================================
import os
import common_functions as cf
	
#===========================================================================================================================
#-------------------------------------------- Simulation Input Results -----------------------------------------------------

#-----------------------------------------------------------------
# Function that prints a list
def print_input_results_list(f,list_name,display_name):
	f.write('\n\n---------------------- '+display_name+' -----------------------\n')
	for name in list_name:
		f.write(str(name)+', ')

#-----------------------------------------------------------------
# Function that prints the MOS Parameters
def print_input_results_mos_parameters(f,optimization_input_parameters):
	f.write('\n\n---------------------- MOS Parameters -----------------------')
	f.write('\nMOSFET File	:'+str(optimization_input_parameters['MOS']['filename']))
	f.write('\nMOS Type 	:'+str(optimization_input_parameters['MOS']['Type']))
	f.write('\nVdd      	:'+str(optimization_input_parameters['MOS']['Vdd']))
	f.write('\nLmin     	:'+str(optimization_input_parameters['MOS']['Lmin']))

#-----------------------------------------------------------------
# Function that prints the output conditions
def print_input_results_output_conditions(f,optimization_input_parameters):
	f.write('\n\n---------------------- Output Conditions -----------------------')
	for name in optimization_input_parameters['output_conditions']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['output_conditions'][name],3))

#-----------------------------------------------------------------
# Function that prints the simulation conditions
def print_input_results_simulation_conditions(f,optimization_input_parameters):
	f.write('\n\n---------------------- Simulation Conditions -----------------------')
	f.write('\nDirectory      :'+str(optimization_input_parameters['simulation']['directory']))
	f.write('\nBasic Filename :'+str(optimization_input_parameters['simulation']['basic_circuit']))
	f.write('\nIIP3 Filename  :'+str(optimization_input_parameters['simulation']['iip3_circuit']))
	f.write('\nTCSH Filename  :'+str(optimization_input_parameters['simulation']['tcsh']))
	f.write('\nStandard Temp  :'+str(optimization_input_parameters['simulation']['std_temp']))
	f.write('\nPin Fixed      :'+str(optimization_input_parameters['simulation']['pin_fixed']))
	f.write('\nPin Start      :'+str(optimization_input_parameters['simulation']['pin_start']))
	f.write('\nPin Stop       :'+str(optimization_input_parameters['simulation']['pin_stop']))
	f.write('\nPin Points     :'+str(optimization_input_parameters['simulation']['pin_points']))
	f.write('\nIIP3 Calculation Points :'+str(optimization_input_parameters['simulation']['iip3_calc_points']))
	
	for name in optimization_input_parameters['simulation']['parameters_list']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['simulation']['parameters_list'][name],3))
	
	for name in optimization_input_parameters['simulation']['cir_writing_dict']:
		f.write('\n'+str(name)+': '+str(optimization_input_parameters['simulation']['cir_writing_dict'][name]))

#-----------------------------------------------------------------
# Function that prints the pre optimization parameters
def print_input_results_pre_optimization(f,optimization_input_parameters):
	f.write('\n\n---------------------- Pre Optimization -----------------------')
	
	f.write('\nStep1b_Limit :'+str(optimization_input_parameters['pre_optimization']['Step1b_Limit']))
	f.write('\nStep2_Limit  :'+str(optimization_input_parameters['pre_optimization']['Step2_Limit']))
	f.write('\nvdsat_reqd      :'+str(optimization_input_parameters['pre_optimization']['vdsat_reqd']))

	f.write('\nPre_Opt_Type	   :'+str(optimization_input_parameters['pre_optimization']['type']))
	f.write('\ngmrs_threshold  :'+str(optimization_input_parameters['pre_optimization']['gmrs_threshold']))
	f.write('\nvdsat_threshold :'+str(optimization_input_parameters['pre_optimization']['vdsat_threshold']))
	
	f.write('\nC1_threshold    :'+str(optimization_input_parameters['pre_optimization']['C1_threshold']))
	f.write('\nC1_threshold    :'+str(optimization_input_parameters['pre_optimization']['C2_threshold']))
	f.write('\nRbias_threshold :'+str(optimization_input_parameters['pre_optimization']['Rbias_threshold']))

	print_input_results_manual_parameters(f,optimization_input_parameters)

#-----------------------------------------------------------------
# Function that prints the manual circuit parameters
def print_input_results_manual_parameters(f,optimization_input_parameters):
	f.write('\nManual Circuit Parameters')
	for name in optimization_input_parameters['pre_optimization']['manual_circuit_parameters']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['pre_optimization']['manual_circuit_parameters'][name],3))

#-----------------------------------------------------------------
# Function that prints the optimization parameters
def print_input_results_optimization(f,optimization_input_parameters):
	f.write('\n\n---------------------- Optimization Parameters -----------------------')
	
	f.write('\nMax Iterations :'+str(optimization_input_parameters['optimization']['max_iteration']))
	f.write('\nAlpha Min      :'+str(optimization_input_parameters['optimization']['alpha_min']))
	f.write('\nConsec Iter    :'+str(optimization_input_parameters['optimization']['consec_iter']))
	
	f.write('\nAlpha Mult      :'+str(optimization_input_parameters['optimization']['alpha_mult']))
	f.write('\nDelta Threshold :'+str(optimization_input_parameters['optimization']['delta_threshold']))
	f.write('\nLoss Type       :'+str(optimization_input_parameters['optimization']['loss_type']))
	f.write('\nUpdate Check    :'+str(optimization_input_parameters['optimization']['update_check']))

	f.write('\nOptimization Name :'+str(optimization_input_parameters['optimization']['optimization_name']))
	f.write('\nOptimization Type :'+str(optimization_input_parameters['optimization']['optimization_type']))

	f.write('\nOptimization Parameters : ')
	for name in optimization_input_parameters['optimization']['optimizing_parameters']:
		f.write(str(name)+' ,')

	print_input_results_loss_weights(f,optimization_input_parameters)
	print_input_results_alpha_parameters(f,optimization_input_parameters)

#-----------------------------------------------------------------
# Function that prints the loss weights
def print_input_results_loss_weights(f,optimization_input_parameters):
	f.write('\n\n---------------------- Loss Weights -----------------------')
	for name in optimization_input_parameters['optimization']['loss_weights']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['optimization']['loss_weights'][name],3))

#-----------------------------------------------------------------
# Function that prints the alpha parameters
def print_input_results_alpha_parameters(f,optimization_input_parameters):
    f.write('\n\n---------------------- Alpha Parameters -----------------------')
    for name in optimization_input_parameters['optimization']['alpha']['values']:
        f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['optimization']['alpha']['values'][name],3))
    f.write('\nAlpha Type  :'+str(optimization_input_parameters['optimization']['alpha']['type']))
    f.write('\nAlpha Start :'+str(optimization_input_parameters['optimization']['alpha']['start']))
    f.write('\nAlpha End   :'+str(optimization_input_parameters['optimization']['alpha']['end']))

#-----------------------------------------------------------------
# Function that prints the acceptable solution
def print_input_acceptable_solution(f,optimization_input_parameters):
	if 'acceptable_solution' not in optimization_input_parameters:
		return
	f.write('\n\n---------------------- Acceptable Solution Parameters -----------------------')
	for name in optimization_input_parameters['acceptable_solution']:
		f.write('\n'+str(name)+': '+cf.num_trunc(optimization_input_parameters['acceptable_solution'][name],3))

#-----------------------------------------------------------------
# Function that prints the filenames
def print_input_results_filenames(f,optimization_input_parameters):
	f.write('\n\n---------------------- Filenames -----------------------')
	f.write('\nRun Status : '+str(optimization_input_parameters['filename']['run_status']))
	f.write('\nOutput     : '+str(optimization_input_parameters['filename']['output']))

#-----------------------------------------------------------------
# Function that prints additional information
def print_input_extra_notes(f):
	f.write('\n\n------------- Note ------------------')
	f.write('\nLoss Type is 1 for normal gradient descent')
	f.write('\nLoss Type is 0 for gradient descent with slope of Io is only considered when other losses are 0 and Io slope is ignored otherwise')
	f.write('\n\nUpdate Check is 1 if we want to perform next iteration with the a previous result having the smaller loss')
	f.write('\nUpdate Check is 0 fif we will perform next iteration with present circuit parameters')

#-----------------------------------------------------------------
# Function that stores input and output data of the simulation
def save_input_results(optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	filename=filename+str('/input_data.txt')
	f=open(filename,'w')

	print_input_results_mos_parameters(f,optimization_input_parameters)
	print_input_results_output_conditions(f,optimization_input_parameters)
	print_input_results_simulation_conditions(f,optimization_input_parameters)
	print_input_results_pre_optimization(f,optimization_input_parameters)
	print_input_results_optimization(f,optimization_input_parameters)
	print_input_acceptable_solution(f,optimization_input_parameters)
	print_input_results_filenames(f,optimization_input_parameters)
	print_input_extra_notes(f)
	
	f.close()
	
#===========================================================================================================================
#-------------------------------------------- Simulation Output Results ----------------------------------------------------

#-----------------------------------------------------------------
# Function that prints parameters
def print_output_parameters(f,parameters):
	for param_name in parameters:
		f.write('\n'+str(param_name)+': '+cf.num_trunc(parameters[param_name],3))

#-----------------------------------------------------------------
# Function that prints parameters
def print_output_mos_parameters(f,optimization_results):
	f.write('\n\n---------------- MOS Parameters ----------------------------------')
	for param_name in optimization_results['mos_parameters']:
		f.write('\n'+str(param_name)+': '+cf.num_trunc(optimization_results['mos_parameters'][param_name],3))

#-----------------------------------------------------------------
# Function that stores output data of the MOS File Calculations
def save_mos_results(mos_parameters,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']
	newpath =filename+'/'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
		
	filename=filename+str('/output_data.txt')
	f=open(filename,'w')

	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOS Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	for param_name in mos_parameters:
		f.write('\n'+str(param_name)+': '+cf.num_trunc(mos_parameters[param_name],3))
	
	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the pre optimization
def save_pre_opt_results(optimization_results,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']
	filename=filename+str('/output_data.txt')
	f=open(filename,'a')

	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	if 'manual_hc' in optimization_results:
		f.write('\n\n--------------------- Manual Hand Calculations ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['manual_hc']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['manual_hc']['extracted_parameters'])

	if 'auto_hc' in optimization_results:
		f.write('\n\n--------------------- Automatic Hand Calculations ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['auto_hc']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['auto_hc']['extracted_parameters'])

	if 'hc_update' in optimization_results:
		f.write('\n\n--------------------- Hand Calculations Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['hc_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['hc_update']['extracted_parameters'])

	if 'gm_update' in optimization_results:
		f.write('\n\n--------------------- gm Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['gm_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['gm_update']['extracted_parameters'])

	if 'gmvd_update' in optimization_results:
		f.write('\n\n--------------------- gmvd Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['gmvd_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['gmvd_update']['extracted_parameters'])
	
	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the optimization
def save_opt_results(optimization_results,optimization_input_parameters):
	filename=optimization_input_parameters['filename']['output']
	filename=filename+str('/output_data.txt')
	f=open(filename,'a')

	print_dict=optimization_results['optimized_results']
	iter_number=print_dict['iter_number']-1

	run_number=optimization_results['run_number']

	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Optimization '+str(run_number)+'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	
	if 'optimization_start' in optimization_results:
		f.write('\n\n--------------------- Optimization Start ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		print_output_parameters(f,optimization_results['optimization_start']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		print_output_parameters(f,optimization_results['optimization_start']['extracted_parameters'])
	

	f.write('\n-------------------------------------------------------------------')
	if optimization_input_parameters['optimization']['optimization_name']=='loss1':
		f.write('\nMaximum Loss of gain+Io+s11+iip3='+cf.num_trunc(print_dict['loss_max'],3))
		f.write('\nOptimized Point occured at iteration='+str(print_dict['iter_number']))
		f.write('\nOptimized Io Loss='+cf.num_trunc(print_dict['Io_loss'],3))
	
	elif optimization_input_parameters['optimization']['optimization_name']=='fom1':
		f.write('\nMaximum Loss of s11='+cf.num_trunc(print_dict['loss_max'],3))
		f.write('\nOptimized Point occured at iteration='+str(print_dict['iter_number']))
		f.write('\nOptimized FOM in dB='+cf.num_trunc(print_dict['FOM'],3))
	
	f.write('\n\n------------------------- Circuit Parameter Values ----------------------------------------')
	print_output_parameters(f,optimization_results['circuit_parameters_iter'][iter_number])
	
	f.write('\n\n------------------------- Output Parameter Values ----------------------------------------')
	print_output_parameters(f,optimization_results['output_parameters_iter'][iter_number])

	if 'acceptable_solution' in optimization_results:
		f.write('Acceptable Solutions:\n')
		for i in optimization_results['acceptable_solution']:
			f.write(str(i)+' ; ')
	
	f.close()


#-----------------------------------------------------------------
# Function that stores output data of the simulation
def save_time_results(timing_results,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	
	filename=filename+str('/output_data.txt')
	f=open(filename,'a')

	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Timing Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	for optimization_name in timing_results:
		f.write('\n\n---------- '+str(optimization_name)+' ----------')
		
		if 'start' in timing_results[optimization_name]:
			f.write('\nStart    : '+str(timing_results[optimization_name]['start']))
			f.write('\nEnd      : '+str(timing_results[optimization_name]['stop']))
			f.write('\nDuration : '+str(timing_results[optimization_name]['stop']-timing_results[optimization_name]['start']))
		
		else:
			for run_number in timing_results[optimization_name]:
				f.write('\n\nRun Number : '+str(run_number))
				f.write('\nStart    : '+str(timing_results[optimization_name][run_number]['start']))
				f.write('\nEnd      : '+str(timing_results[optimization_name][run_number]['stop']))
				f.write('\nDuration : '+str(timing_results[optimization_name][run_number]['stop']-timing_results[optimization_name][run_number]['start']))
	
	f.close()


#===========================================================================================================================

