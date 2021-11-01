#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: pre_optimization.py
File Description 	: This file will perform pre optimization and calculate an initial point to be used 
					  at the start of the gradient descent algorithm

Functions structure in this file:
	--> pre_optimization
		--> save_input_results_pre_optimization
		--> save_mos_results
		--> manual_circuit_parameters
		--> calculate_mos_parameters
		--> save_output_results_pre_optimization

"""

#===========================================================================================================================
import datetime
import CG_LNA.extra_function as cff
#import CG_LNA.spectre as sp
import CG_LNA.hand_calculation_1 as hc1
import CG_LNA.hand_calculation_2 as hc2


#===========================================================================================================================
#----------------------------------- File Writing Functions ----------------------------------------------------------------

#-----------------------------------------------------------------
# Function that stores input data of the simulation
# Inputs  : optimization_input_parameters
# Outputs : NONE
def save_input_results_pre_optimization(optimization_input_parameters):

	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/input_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	f.write('\nStep1b_Limit :'+str(optimization_input_parameters['pre_optimization']['Step1b_Limit']))
	f.write('\nStep2_Limit  :'+str(optimization_input_parameters['pre_optimization']['Step2_Limit']))
	f.write('\nvdsat_reqd      :'+str(optimization_input_parameters['pre_optimization']['vdsat_reqd']))

	f.write('\nPre_Opt_Type	   :'+str(optimization_input_parameters['pre_optimization']['type']))
	f.write('\ngmrs_threshold  :'+str(optimization_input_parameters['pre_optimization']['gmrs_threshold']))
	f.write('\nvdsat_threshold :'+str(optimization_input_parameters['pre_optimization']['vdsat_threshold']))
	
	f.write('\nC1_threshold    :'+str(optimization_input_parameters['pre_optimization']['C1_threshold']))
	f.write('\nC1_threshold    :'+str(optimization_input_parameters['pre_optimization']['C2_threshold']))
	f.write('\nRbias_threshold :'+str(optimization_input_parameters['pre_optimization']['Rbias_threshold']))

	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the pre optimization
# Inputs  : optimization_results, optimization_input_parameters
# Outputs : NONE
def save_output_results_pre_optimization(optimization_results,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	
	if 'manual_hc' in optimization_results:
		f.write('\n\n--------------------- Manual Hand Calculations ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['manual_hc']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['manual_hc']['extracted_parameters'])

	if 'auto_hc' in optimization_results:
		f.write('\n\n--------------------- Automatic Hand Calculations ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['auto_hc']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['auto_hc']['extracted_parameters'])

	if 'hc_update' in optimization_results:
		f.write('\n\n--------------------- Hand Calculations Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['hc_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['hc_update']['extracted_parameters'])

	if 'gm_update' in optimization_results:
		f.write('\n\n--------------------- gm Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['gm_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['gm_update']['extracted_parameters'])

	if 'gmvd_update' in optimization_results:
		f.write('\n\n--------------------- gmvd Update ---------------------------------')
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['gmvd_update']['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cff.print_output_parameters(f,optimization_results['gmvd_update']['extracted_parameters'])
	
	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the MOS File Calculations
# Inputs  : mos_parameters, optimization_input_parameters
# Outputs : NONE
def save_mos_results(mos_parameters,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ MOS Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	for param_name in mos_parameters:
		f.write('\n'+str(param_name)+': '+cff.num_trunc(mos_parameters[param_name],3))
	
	f.close()


#===========================================================================================================================
#----------------------------------- Defining the functions for simple calculations ----------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Function to manually choose the Initial Circuit Parameters
# Inputs  : optimization_input_parameters
# Outputs :	circuit_parameters, extracted_parameters
def manual_initial_parameters(cir,optimization_input_parameters):

	# Getting Circuit Parameters
	circuit_parameters=optimization_input_parameters['pre_optimization']['manual_circuit_parameters'].copy()
		
	# Running Eldo
	#extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)	

	extracted_parameters=cir.update_circuit(circuit_parameters)
	
	return circuit_parameters,extracted_parameters


#===========================================================================================================================
#------------------------------------ MOSFET EXTRACTION --------------------------------------------------------------------

#-----------------------------------------------------------------
# Function that extracts the MOSFET File Parameeters
# Inputs  : Optimization Input Parameters
# Outputs : MOS_Parameters
def calculate_mos_parameters(optimization_input_parameters):
	
	# Setting Lmin and Vdd
	Lmin=optimization_input_parameters['MOS']['Lmin']
	vdd=optimization_input_parameters['MOS']['Vdd']
	cox=optimization_input_parameters['MOS']['cox']
	un=optimization_input_parameters['MOS']['un']
	vt=optimization_input_parameters['MOS']['vt']

	# Extracting From File
	mos_parameters = {'un':un,'cox':cox,'vt':vt,'Lmin':Lmin,'vdd':vdd}
	
	# Printing the MOSFET Parameters
	cff.print_MOS_parameters(mos_parameters)

	# Storing the results
	save_mos_results(mos_parameters,optimization_input_parameters)

	return mos_parameters


#===========================================================================================================================
#------------------------------------------- Output Functions --------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Function to perform pre-optimization
# Inputs  : mos_parameters, optimization_input_parameters, timing_results
# Outputs :	circuit_parameters, extracted_parameters
def pre_optimization(cir,optimization_input_parameters,timing_results):

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Pre Optimization Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['pre_optimization']={}
	timing_results['pre_optimization']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('*********************************** Pre Optimization *******************************************************')

	save_input_results_pre_optimization(optimization_input_parameters)

	cff.write_simulation_parameters(optimization_input_parameters,'pre_optimization',0)
	cir.optimization_input_parameters=optimization_input_parameters

	optimization_results={}
	
	#======================================================== Manual Initial Points =============================================================================================================

	if optimization_input_parameters['pre_optimization']['type']=='manual':
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manual Operating Point Selection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#--------------------Initial Point Calculations-------------------------

		# Calculating the Values of Circuit Parameters
		circuit_parameters,extracted_parameters=manual_initial_parameters(cir,optimization_input_parameters)

		# Storing the Circuit and Extracted Parameters
		optimization_results['manual_hc']={}
		optimization_results['manual_hc']['circuit_parameters']=circuit_parameters.copy()
		optimization_results['manual_hc']['extracted_parameters']=extracted_parameters.copy()

		# Printing the values
		cff.print_circuit_parameters(circuit_parameters)
		cff.print_extracted_outputs(extracted_parameters)

		#cff.wait_key()

	
	#======================================================== Automatic Initial Points =============================================================================================================

	if optimization_input_parameters['pre_optimization']['type']==1:
		
		print('************************************************************************************************************')
		print('*********************************** MOSFET Parameters ******************************************************')

		# Extracting the MOSFET Parameters from the MOS file
		mos_parameters=calculate_mos_parameters(optimization_input_parameters)
		
		circuit_parameters,extracted_parameters=hc1.automatic_initial_parameters(cir,mos_parameters,optimization_input_parameters,optimization_results)

	if optimization_input_parameters['pre_optimization']['type']==2:

		print('************************************************************************************************************')
		print('*********************************** MOSFET Parameters ******************************************************')

		# Extracting the MOSFET Parameters from the MOS file
		mos_parameters=calculate_mos_parameters(optimization_input_parameters)

		circuit_parameters,extracted_parameters=hc2.automatic_initial_parameters(cir,mos_parameters,optimization_input_parameters,optimization_results)

	# Printing the values
	cff.print_circuit_parameters(circuit_parameters)
	cff.print_extracted_outputs(extracted_parameters)

	# Storing the results
	save_output_results_pre_optimization(optimization_results,optimization_input_parameters)

	# Storing the finishing time
	timing_results['pre_optimization']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Pre Optimization End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	return circuit_parameters,extracted_parameters
	

#===========================================================================================================================
