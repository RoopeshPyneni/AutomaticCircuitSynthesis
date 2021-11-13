#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: CS_LNA/hand_calculation_2.py
File Description 	: This file will initialize the optimization_input_parameters and run the complete_optimization file

Functions structure in this file:
	--> TBD


"""

#===========================================================================================================================
import numpy as np
import CS_LNA.extra_function as cff # type: ignore


"""
===========================================================================================================================
------------------------------------Defining the functions for simple calculations-----------------------------------------
"""


#-----------------------------------------------------------------------------------------------
# Converting from dB to normal
# Outputs : normal values
def db_to_normal(val_db):
    return 10**(val_db/10)

#-----------------------------------------------------------------------------------------------
# Calculating the resistance from the inductance values
# Outputs : resistance
def calculate_resistance_inductor(ind,fo,Q):
    res=Q*2*np.pi*fo*ind
    return res

#-----------------------------------------------------------------------------------------------
# Calculating Ld and Rd
# Outputs : Ld, Rd
def calculate_Ld(Cload,fo):
    wo=2*np.pi*fo
    Ld=1/(wo*wo*Cload)
    return Ld

#-----------------------------------------------------------------------------------------------
# Calculating gm
# Outputs : gm
def calculate_gm(gain,Rd):
    return gain/(2*Rd*2.5)

#-----------------------------------------------------------------------------------------------
# Calculating Cgs
# Outputs : cgs
def calculate_cgs(Rs,fo):
    wo=2*np.pi*fo
    cgs=1/(2*Rs*2.5*wo)
    return cgs

#-----------------------------------------------------------------------------------------------
# Calculating W
# Outputs : W
def calculate_W(cgs,Lmin,Cox):
    return 3*cgs/(2*Lmin*Cox)

#-----------------------------------------------------------------------------------------------
# Calculating C1
# Outputs : C1
def calculate_C1(Cox,Lmin,W,cgs):
    return cgs-(2/3*Cox*Lmin*W)

#-----------------------------------------------------------------------------------------------
# Calculating Io
# Outputs : Io
def calculate_Io(gm,un,cox,W,Lmin):
    return (gm*gm)/(2*un*cox*W/Lmin)

#-----------------------------------------------------------------------------------------------
# Calculating Ls and Rs
# Outputs : Ls,Rs
def calculate_Ls(rsource,cgs,gm,fo):
    return rsource*cgs/gm
    
#-----------------------------------------------------------------------------------------------
# Calculating Lg and Rg
# Outputs : Lg,Rg
def calculate_Lg(Ls,cgs,fo):
    w=2*np.pi*fo
    Lg=1/(w*w*cgs)-Ls
    return Lg
	

	

"""
===========================================================================================================================
-------------------------------------------- Main Functions ---------------------------------------------------------------
"""

	
#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the Initial Circuit Parameters	
# Inputs  : mos_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, extracted_parameters
def calculate_initial_parameters(cir,optimization_input_parameters):

	output_conditions=optimization_input_parameters['output_conditions']
    
    # Getting the output conditions
	Cload=output_conditions['Cload']
	fo=output_conditions['wo']/(2*np.pi)
	Rs=output_conditions['Rs']
	Lmin=cir.mos_parameters['Lmin']
	Cox=cir.mos_parameters['cox']
	un=cir.mos_parameters['un']

    # Calculating the circuit parameters
	circuit_parameters={}
	circuit_parameters['Ld']=calculate_Ld(Cload,fo)
	gm=20e-3
	cgs=calculate_cgs(Rs,fo)
	circuit_parameters['W']=calculate_W(cgs,Lmin,Cox)
	circuit_parameters['Io']=calculate_Io(gm,un,Cox,circuit_parameters['W'],Lmin)
	circuit_parameters['Ls']=calculate_Ls(Rs,cgs,gm,fo)
	circuit_parameters['Lg']=calculate_Lg(circuit_parameters['Ls'],cgs,fo)
	circuit_parameters['Rb']=5000
	circuit_parameters['Cs']=100/(2*np.pi*50*fo)

	# Running the circuit
	cir.update_circuit(circuit_parameters)


#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters	after calculating the new value of vt
# Inputs  : circuit_parameters, mos_parameters, extracted_parameters, optimization_input_parameters
# Outputs : circuit_parameters, dc_outputs, mos_parameters, extracted_parameters
def update_initial_parameters(cir,optimization_input_parameters):

	i=0
	while i<5 and cir.extracted_parameters['s11_db']>-15.0:

		# Printing the iteration number
		i+=1
		print('----- Iteration ',i,' -----')

		# Updating the values
		fo=optimization_input_parameters['output_conditions']['wo']/(2*np.pi)
		cir.circuit_parameters['Ls']=cir.circuit_parameters['Ls']*50/cir.extracted_parameters['1_Zin_R']
		cir.circuit_parameters['Lg']=cir.circuit_parameters['Lg']-cir.extracted_parameters['1_Zin_I']/(2*np.pi*fo)
		
		# Running the circuit
		cir.run_circuit()
	


"""
===========================================================================================================================
-------------------------------------------- Output Functions -------------------------------------------------------------
"""


#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the initial parameters by completing all the sub steps of pre optimization
# Inputs  : mos_parameters, optimization_input_parameters, optimization_results
# Outputs : circuit_parameters, extracted_parameters
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
		
	print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


	#======================================================== Step 1 =======================================================
	print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	calculate_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['auto_hc']={}
	optimization_results['auto_hc']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['auto_hc']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cff.print_circuit_parameters(cir.circuit_parameters)
	cff.print_extracted_outputs(cir.extracted_parameters)

	

	#======================================================== Step 2 =======================================================
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	update_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['hc_update']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cff.print_circuit_parameters(cir.circuit_parameters)
	cff.print_extracted_outputs(cir.extracted_parameters)

#===========================================================================================================================
