#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform hand calculations using method 2
					  This method involves choosing the points such that gm is low and Rb is 50
"""

#===========================================================================================================================
import numpy as np
import CG_LNA.extra_function as cff # type: ignore

#===========================================================================================================================
#------------------------------------Defining the functions for simple calculations-----------------------------------------

#-----------------------------------------------------------------------------------------------
# Calculating C1
# Inputs  : output_conditions, optimization_input_parameters
# Outputs : C1
def calc_C1(opt_conditions,optimization_input_parameters):
	
	threshold1=optimization_input_parameters['pre_optimization']['C1_threshold']

	# Assigning the values
	wo=opt_conditions['wo']
	Rs=opt_conditions['Rs']
	
	# Calculating C1
	C1=threshold1/(wo*Rs)
	
	return C1

#-----------------------------------------------------------------------------------------------
# Calculating W
# Inputs  : circuit_parameters, mos_parameters
# Outputs : W
def calc_W(circuit_parameters,mos_parameters):

	# Assigning the values
	gm=2e-3
	Io=circuit_parameters['Io']
	Lmin=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	un=mos_parameters['un']
	
	# Calculating W
	W=(Lmin*gm*gm)/(un*cox*2*Io)
	
	return W

#-----------------------------------------------------------------------------------------------	
# Calculating Rd by assuming gain=Rd/2*Rs
# Inputs  : output_conditions
# Outputs : Rd
def calc_Rd(opt_conditions):

	# Calculating the required gain
	gain=np.sqrt(cff.db_to_normal(opt_conditions['gain_db']))
	gm=2e-3

	# Calculating Rd
	Rd=gain/gm
	
	return Rd

#-----------------------------------------------------------------------------------------------
# Calculating C2 from W and Rbias 
# Inputs  : mos_parameters, output_conditions, circuit_parameters, optimization_input_parameters
# Outputs : C2, Rbias
def calc_C2(mos_parameters,opt_conditions,circuit_parameters,optimization_input_parameters):
	
	threshold2=optimization_input_parameters['pre_optimization']['C2_threshold']
	threshold3=optimization_input_parameters['pre_optimization']['Rbias_threshold']
	Rbias_min=optimization_input_parameters['pre_optimization']['Rbias_minimum']

	# Assigning the values
	L=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	wo=opt_conditions['wo']
	W=circuit_parameters['W']
	
	# Calculating C2
	C2=threshold2*(2*W*L*cox/3)
	
	# Calculating Rbias
	Rbias=threshold3/(wo*C2)
	Rbias_min=max(Rbias,Rbias_min)
	
	return C2,Rbias

#-----------------------------------------------------------------------------------------------
# Calculating DC Output Values from Initial Circuit Parameters
# Inputs  : circuit_parameters,mos_parameters,output_conditions
# Outputs : dc_outputs
def calc_dc_opt(circuit_parameters,mos_parameters,opt_conditions):
	vs=circuit_parameters['Io']*circuit_parameters['Rb']
	vgs=mos_parameters['vt']+ 2*opt_conditions['Rs']*circuit_parameters['Io']
	vg=vgs+vs
	vd=mos_parameters['vdd']-circuit_parameters['Io']*circuit_parameters['Rd']
	
	dc_outputs={'vg':vg,'vd':vd,'vs':vs}
	
	return dc_outputs

#-----------------------------------------------------------------------------------------------
# Calculating C2 from Cgs & Cgd and Rbias
# Inputs  : extracted_parameters, output_conditions, circuit_parameters, optimization_input_parameters
# Outputs : C2, Rbias
def calc_C2_updated(extracted_parameters,opt_conditions,circuit_parameters,optimization_input_parameters):
	
	threshold2=optimization_input_parameters['pre_optimization']['C2_threshold']
	threshold3=optimization_input_parameters['pre_optimization']['Rbias_threshold']
	Rbias_min=optimization_input_parameters['pre_optimization']['Rbias_minimum']

	# Assigning the values
	cgs=extracted_parameters['cgs1']
	cgd=extracted_parameters['cgd1']
	gain=circuit_parameters['Rd']/(2*opt_conditions['Rs'])
	wo=opt_conditions['wo']
	
	# Calculating C2
	C2a=threshold2*cgs
	C2b=threshold2*cgd*gain
	C2=np.maximum(C2a,C2b)
	
	# Calculating Rbias
	Rbias=threshold3/(wo*C2)
	Rbias=max(Rbias,Rbias_min)
	
	return C2,Rbias


	
#===========================================================================================================================
#-------------------------------------------- Main Functions ---------------------------------------------------------------
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the Initial Circuit Parameters
# Inputs  : mos_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, extracted_parameters
def calculate_initial_parameters(cir,optimization_input_parameters):

	opt_conditions=optimization_input_parameters['output_conditions']
	
	# Creating dictionary of output parameters	
	circuit_parameters={}
	circuit_parameters_list=['Rb','Rd','Io','C1','C2','W','Rbias']
	for param_name in circuit_parameters_list:
		circuit_parameters[param_name]=0

	# Calculating Rb
	circuit_parameters['Rb']=optimization_input_parameters['output_conditions']['Rs']

	# Calculating C1
	circuit_parameters['C1']=calc_C1(opt_conditions,optimization_input_parameters)
	
	# Calculating Io
	circuit_parameters['Io']=100e-6

	# Calculating W
	circuit_parameters['W']=calc_W(circuit_parameters,cir.mos_parameters)
	
	# Calculating Rd
	circuit_parameters['Rd']=calc_Rd(opt_conditions)
	
	# Calculating C2 and Rbias
	circuit_parameters['C2'],circuit_parameters['Rbias']=calc_C2(cir.mos_parameters,opt_conditions,circuit_parameters,optimization_input_parameters)
	
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(circuit_parameters,cir.mos_parameters,opt_conditions)
	
	# Running Eldo
	cir.update_circuit(circuit_parameters)
	
	return dc_outputs


#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters	
# Inputs  : circuit_parameters, mos_parameters, extracted_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, mos_parameters, extracted_parameters
def update_initial_parameters(cir,optimization_input_parameters):

	opt_conditions=optimization_input_parameters['output_conditions']
	cir.mos_parameters['vt']=cir.extracted_parameters['vt']
		
	# Calculating C2 and Rbias
	cir.circuit_parameters['C2'],cir.circuit_parameters['Rbias']=calc_C2_updated(cir.extracted_parameters,opt_conditions,cir.circuit_parameters,optimization_input_parameters)
		
	# Running Eldo
	cir.run_circuit()
		
	# Updating the value of vt
	cir.mos_parameters['vt']=cir.extracted_parameters['vt']	
		
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(cir.circuit_parameters,cir.mos_parameters,opt_conditions)
	
	return dc_outputs



#===========================================================================================================================
#-------------------------------------------- Output Functions -------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the initial parameters by completing all the sub steps of pre optimization
# Inputs  : mos_parameters, optimization_input_parameters, optimization_results
# Outputs : circuit_parameters, extracted_parameters
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
	
	print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	#======================================================== Step 1 =============================================================================================================
	print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	dc_initial_outputs=calculate_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['auto_hc']={}
	optimization_results['auto_hc']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['auto_hc']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cff.print_circuit_parameters(cir.circuit_parameters)
	cff.print_DC_outputs(dc_initial_outputs,cir.mos_parameters)
	cff.print_extracted_outputs(cir.extracted_parameters)



	#======================================================== Step 1 b ============================================================================================================
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	dc_initial_outputs=update_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['hc_update']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cff.print_circuit_parameters(cir.circuit_parameters)
	cff.print_DC_outputs(dc_initial_outputs,cir.mos_parameters)
	cff.print_extracted_outputs(cir.extracted_parameters)
	
#===========================================================================================================================
