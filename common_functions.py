#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file contains many commonly used functions
"""

#===========================================================================================================================
import numpy as np
import os
#===========================================================================================================================


"""
===========================================================================================================================
---------------------------------------- Calculation Functions ------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# These functions change the parameters from normal scale to dB and dBm scale and vice versa
def db_to_normal(x_db):
	exp_x=x_db/10.0
	x_normal=10.0**exp_x
	return x_normal

def dbm_to_normal(x_dbm):
	x_db=x_dbm-30.0
	exp_x=x_db/10.0
	x_normal=10.0**exp_x
	return x_normal
	
def normal_to_db(x_normal):
	x_db=10.0*np.log10(x_normal)
	return x_db
	
def normal_to_dbm(x_normal):
	x_db=10.0*np.log10(x_normal)
	x_dbm=x_db+30.0
	return x_dbm

#-----------------------------------------------------------------------------------------------
# This function truncates a number to a few significant numbers and returns a string
def num_trunc(num,pts):
	
	# Initializing Length of number
	num_len=-1 
	
	# Checking for zero
	if num==0:
		return '0'

	# Converting to a positive number
	flag=1
	if num<0:
		num*=-1
		flag=-1
		
	# Truncating the number
	i=-5
	while i<20:
		if 10**(pts-1)<=num*10**i and 10**pts>num*10**i:
			val=int(num*10**i)
			num=val*10**(-1*i)
			break
		else:
			i+=1
	
	# Finding the number and string
	flag_check=1
	
	if num>1:
		return str(flag*num)
	
	dict_char={0:' ',1:'m',2:'u',3:'n',4:'p',5:'f',6:'a'}
	if num>=1:
		val=num
		char=dict_char[0]
		flag_check=0
		
	if flag_check==1:
		i=0
		while i<6:
			num*=1000
			i+=1
			if num>=1 and num<1000:
				val=num
				char=dict_char[i]
				flag_check=0
				break
				
	if val>=100 and val<1000:
		num_len=pts #
	else:
		num_len=(1+pts) #
		
		
	if flag_check==1:
		val=num*1000
		char=dict_char[6]
	
	
	# Converting the number into a complete string
	val=val*flag
	char_val=str(val)
	
	if flag==-1:
		num_len+=1
		
	if num_len>0:
		char_val=char_val[:num_len]
	
	if char!=' ':
		char_val=char_val+' '+char
		
	return char_val


"""
===========================================================================================================================
------------------------------------ Dictionary Modification Functions ----------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function that is used to modify the simulation parameters
def write_simulation_parameters(optimization_input_parameters,optimization_name,iteration_number):
	if iteration_number==0:
		if 'parameters_list' in optimization_input_parameters[optimization_name]['simulation']:
			for param_name in optimization_input_parameters[optimization_name]['simulation']['parameters_list']:
				optimization_input_parameters['simulation']['parameters_list'][param_name]=optimization_input_parameters[optimization_name]['simulation']['parameters_list'][param_name]
		
		for param_name in optimization_input_parameters[optimization_name]['simulation']:
			if param_name != 'parameters_list':
				optimization_input_parameters['simulation'][param_name]=optimization_input_parameters[optimization_name]['simulation'][param_name]
	
	else:
		if 'parameters_list' in optimization_input_parameters[optimization_name]['simulation'][iteration_number]:
			for param_name in optimization_input_parameters[optimization_name]['simulation'][iteration_number]['parameters_list']:
				optimization_input_parameters['simulation']['parameters_list'][param_name]=optimization_input_parameters[optimization_name]['simulation'][iteration_number]['parameters_list'][param_name]
		
		for param_name in optimization_input_parameters[optimization_name]['simulation'][iteration_number]:
			if param_name != 'parameters_list':
				optimization_input_parameters['simulation'][param_name]=optimization_input_parameters[optimization_name]['simulation'][iteration_number][param_name]

#-----------------------------------------------------------------------------------------------
# This function is used to update the MOS Parameters
def update_MOS_parameters(mos,un,cox,vt,Lmin,vdd):
	mos['un']=un
	mos['cox']=cox
	mos['vt']=vt
	mos['un']=Lmin
	mos['un']=vdd
	return mos


"""
===========================================================================================================================
--------------------------------------- Output Printing Functions ---------------------------------------------------------
"""

trunc_val=3

#-----------------------------------------------------------------
# Function that prints parameters completely without truncation
def print_output_parameters_complete(f,parameters):
	for param_name in parameters:
		f.write('\n'+str(param_name)+': '+str(parameters[param_name]))
		
#-----------------------------------------------------------------
# Function that prints parameters
def print_output_parameters(f,parameters):
	for param_name in parameters:
		f.write('\n'+str(param_name)+': '+num_trunc(parameters[param_name],3))

#-----------------------------------------------------------------
# Function that prints simulation parameters
def print_simulation_parameters(f,cir):

	f.write('\n\n')
	f.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	for param_name in cir.circuit_initialization_parameters['simulation']['standard_parameters']:
		f.write('\n'+str(param_name)+': '+str(cir.circuit_initialization_parameters['simulation']['standard_parameters'][param_name]))

#-----------------------------------------------------------------------------------------------
# This function is used to wait for key press
def wait_key():
	input('\n\nPress Enter to continue')
	os.system("clear")
	
#-----------------------------------------------------------------------------------------------
# Printing the MOSFET Parameters
def print_MOS_parameters(mos_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------------MOSFET Parameters--------------------------\n')
	print ('uncox = ', num_trunc(mos_parameters['un']*mos_parameters['cox'],trunc_val))
	print ('un    = ', num_trunc(mos_parameters['un'],trunc_val))
	print ('cox   = ', num_trunc(mos_parameters['cox'],trunc_val))
	print ('vt    = ', num_trunc(mos_parameters['vt'],trunc_val))
	print ('Lmin  = ', num_trunc(mos_parameters['Lmin'],trunc_val))
	print ('vdd   = ', num_trunc(mos_parameters['vdd'],trunc_val))
	
#-----------------------------------------------------------------------------------------------
# Printing the circuit parameters
def print_circuit_parameters(circuit_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------------Circuit Parameters-------------------------\n')
	for param_name in circuit_parameters:
		print(param_name,' = ',num_trunc(circuit_parameters[param_name],trunc_val))
		
#-----------------------------------------------------------------------------------------------
# Printing the extracted parameters
def print_extracted_parameters(extracted_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------------Extracted Parameters-----------------------\n')
	for param_name in extracted_parameters:
		print(param_name,' = ',num_trunc(extracted_parameters[param_name],trunc_val))
	
#-----------------------------------------------------------------------------------------------
# Printing the loss parameters
def print_loss_parameters(loss_parameters):
	for param_name in loss_parameters:
		print(param_name,' = ',num_trunc(loss_parameters[param_name],trunc_val))
	
#===========================================================================================================================
