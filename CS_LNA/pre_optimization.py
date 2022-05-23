#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform pre optimization and calculate an initial point to be used 
					  at the start of the gradient descent algorithm
"""

#===========================================================================================================================
import datetime
import common_functions as cf		# type: ignore
import CS_LNA.hand_calculation_1 as hc1 # type: ignore
import CS_LNA.hand_calculation_2 as hc2 # type: ignore
import CS_LNA.hand_calculation_3 as hc3 # type: ignore
import CS_LNA.hand_calculation_4 as hc4 # type: ignore

"""
===========================================================================================================================
----------------------------------- FILE WRITING FUNCTIONS ----------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that stores input data of the simulation
def save_input_results_pre_optimization(cir,optimization_input_parameters):

	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/input_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n')
	f.write('********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization Parameters ~~~~~~~~~~~~~~~~~~~~~~~~\n')
	
	f.write('\nPre_Opt_Type		:'+str(optimization_input_parameters['pre_optimization']['type']))
	f.write('\nI_max R_divider	:'+str(optimization_input_parameters['pre_optimization']['I_Rdivider_max']))

	cf.print_simulation_parameters(f,cir)
	
	f.close()

#-----------------------------------------------------------------
# Function that stores output data of the pre optimization
def save_output_results_pre_optimization(optimization_results,optimization_input_parameters):
	
	# Opening the file
	filename=optimization_input_parameters['filename']['output']+str('/output_data.txt')
	f=open(filename,'a')

	# Storing the results
	f.write('\n\n')
	f.write('********************************************************************************\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pre Optimization Parameters ~~~~~~~~~~~~~~~~~~~~~~~~\n')
	
	for name in optimization_results:
		f.write('\n\n--------------------- '+str(name)+' ---------------------------------')
		f.write('\n\n---------------- Initial Circuit Parameters ---------------')
		cf.print_output_parameters_complete(f,optimization_results[name]['initial_circuit_parameters'])
		f.write('\n\n---------------- Circuit Parameters ------------------------')
		cf.print_output_parameters(f,optimization_results[name]['circuit_parameters'])
		f.write('\n\n---------------- Extracted Parameters ------------------------')
		cf.print_output_parameters(f,optimization_results[name]['extracted_parameters'])

	f.close()


"""
===========================================================================================================================
----------------------------------- MANUAL HAND CALCULATIONS --------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to manually choose the Initial Circuit Parameters
def manual_initial_parameters(cir,optimization_input_parameters):

	# Running Eldo
	cir.update_circuit(optimization_input_parameters['pre_optimization']['manual_circuit_parameters'].copy())


"""
===========================================================================================================================
----------------------------------- OUTPUT FUNCTIONS ----------------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------
# Function to perform pre-optimization
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
	print('')
	
	cir.update_simulation_parameters(optimization_input_parameters['pre_optimization']['simulation'])
	save_input_results_pre_optimization(cir,optimization_input_parameters)

	optimization_results={}
	
	#=============================== Manual Initial Points ================================================================

	if optimization_input_parameters['pre_optimization']['type']=='manual':
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Manual Operating Point Selection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#--------------------Initial Point Calculations-------------------------

		# Calculating the Values of Circuit Parameters
		manual_initial_parameters(cir,optimization_input_parameters)

		# Storing the Circuit and Extracted Parameters
		optimization_results['manual_hc']={}
		optimization_results['manual_hc']['initial_circuit_parameters']=cir.get_initial_circuit_parameters()
		optimization_results['manual_hc']['circuit_parameters']=cir.get_circuit_parameters()
		optimization_results['manual_hc']['extracted_parameters']=cir.get_extracted_parameters()

	#=============================== Automatic Initial Points =============================================================
	
	if optimization_input_parameters['pre_optimization']['type']==1:
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')

		# Extracting the MOSFET Parameters from the MOS file
		hc1.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)
	
	if optimization_input_parameters['pre_optimization']['type']==2:
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')

		# Extracting the MOSFET Parameters from the MOS file
		hc2.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)
	
	if optimization_input_parameters['pre_optimization']['type']==3:
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')

		# Extracting the MOSFET Parameters from the MOS file
		hc3.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)

	if optimization_input_parameters['pre_optimization']['type']==4:
		
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')

		# Extracting the MOSFET Parameters from the MOS file
		hc4.automatic_initial_parameters(cir,optimization_input_parameters,optimization_results)

	# Printing the values
	cf.print_initial_circuit_parameters(cir.get_initial_circuit_parameters())
	cf.print_circuit_parameters(cir.get_circuit_parameters())
	cf.print_extracted_parameters(cir.get_extracted_parameters())

	# Storing the results
	save_output_results_pre_optimization(optimization_results,optimization_input_parameters)

	# Storing the finishing time
	timing_results['pre_optimization']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Pre Optimization End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	

#===========================================================================================================================
