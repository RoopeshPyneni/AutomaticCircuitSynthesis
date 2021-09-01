#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Common Functions File:
"""
#===========================================================================================================================
import numpy as np
import os
#===========================================================================================================================



#===========================================================================================================================
#---------------------------------------- Calculation Functions ------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#These functions change the parameters from normal scale to dB and dBm scale and vice versa
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
# This function returns the truncated number as a string
def num_trunc(num,pts):
	
	num_len=-1 #
	
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


#===========================================================================================================================
#------------------------------------ Dictionary Modification Functions ----------------------------------------------------

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

	
#===========================================================================================================================
#--------------------------------------- Output Printing Functions ---------------------------------------------------------

trunc_val=3

#-----------------------------------------------------------------------------------------------
def wait_key():
	input('\n\nPress Enter to continue')
	os.system("clear")


#-----------------------------------------------------------------------------------------------
#Assigning MOSFET Parameters
def update_MOS_parameters(mos,un,cox,vt,Lmin,vdd):
	mos['un']=un
	mos['cox']=cox
	mos['vt']=vt
	mos['un']=Lmin
	mos['un']=vdd
	return mos
	
"""
#-----------------------------------------------------------------------------------------------
# Printing the MOSFET Parameters
def print_change_opt_gm(W1,Io1,Rb1,W2,Io2,Rb2):
	print ('\n____________________________________________________________________')
	print('------------------------Changes in parameters------------------------\n')
	print('Previous Values:')
	print('W   = ',num_trunc(W1 ,3))
	print('Io  = ',num_trunc(Io1,3))
	print('Rb  = ',num_trunc(Rb1,3))
	print('\nUpdated Values:')
	print('W   = ',num_trunc(W2 ,3))
	print('Io  = ',num_trunc(Io2,3))
	print('Rb  = ',num_trunc(Rb2,3))
"""


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
	print ('W     = ', num_trunc(circuit_parameters['W'],trunc_val))
	print ('Rb    = ',num_trunc(circuit_parameters['Rb'],trunc_val))
	print ('Rd    = ',num_trunc(circuit_parameters['Rd'],trunc_val))
	print ('Io    = ',num_trunc(circuit_parameters['Io'],trunc_val))
	print ('C1    = ',num_trunc(circuit_parameters['C1'],trunc_val))
	print ('C2    = ',num_trunc(circuit_parameters['C2'],trunc_val))
	print ('Rbias = ', num_trunc(circuit_parameters['Rbias'],trunc_val))
	
#-----------------------------------------------------------------------------------------------
# Printing the DC outputs
def print_DC_outputs(dc_outputs,mos_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------Hand Calculation DC Outputs----------------------\n')
	print ('vg     = ',num_trunc(dc_outputs['vg'],trunc_val))
	print ('vs     = ',num_trunc(dc_outputs['vs'],trunc_val))
	print ('vd     = ',num_trunc(dc_outputs['vd'],trunc_val))
	print ('vgs-vt = ',num_trunc(dc_outputs['vg']-dc_outputs['vs']-mos_parameters['vt'],trunc_val))
		
#-----------------------------------------------------------------------------------------------
# Printing the extracted parameters
def print_extracted_outputs(extracted_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------------Extracted Outputs--------------------------\n')
	print ('gm1      = ',num_trunc(extracted_parameters['gm1'],trunc_val))
	print ('gds1     = ',num_trunc(extracted_parameters['gds1'],trunc_val))
	print ('vt1      = ',num_trunc(extracted_parameters['vt'],trunc_val))
	
	print ('vd       = ',num_trunc(extracted_parameters['vd'],trunc_val))
	print ('vg       = ',num_trunc(extracted_parameters['vg'],trunc_val))
	print ('vs       = ',num_trunc(extracted_parameters['vs'],trunc_val))
	print ('Io       = ',num_trunc(extracted_parameters['Io'],trunc_val))
	
	print ('2*Io/gm  = ',num_trunc(2*extracted_parameters['Io']/extracted_parameters['gm1'],trunc_val))
	print ('vgs-vt   = ',num_trunc(extracted_parameters['vg']-extracted_parameters['vs']-extracted_parameters['vt'],trunc_val))
	print ('vdsat    = ',num_trunc(extracted_parameters['vdsat'],trunc_val))
	
	print ('cgs1     = ',num_trunc(extracted_parameters['cgs1'],trunc_val))
	print ('cgd1     = ',num_trunc(extracted_parameters['cgd1'],trunc_val))
	
	print ('iip3_dbm = ',num_trunc(extracted_parameters['iip3_dbm'],trunc_val))
	print ('s11_db   = ',num_trunc(extracted_parameters['s11_db'],trunc_val))
	print ('Gain_db  = ',num_trunc(extracted_parameters['gain_db'],trunc_val))
	print ('NF_db    = ',num_trunc(extracted_parameters['nf_db'],trunc_val))
	
	
#-----------------------------------------------------------------------------------------------
# Printing the extracted parameters for main optimization
def print_extracted_outputs_optimization(extracted_parameters):
	print ('\n____________________________________________________________________')
	print ('-------------------------Extracted Outputs--------------------------\n')
	print ('Io       = ',num_trunc(extracted_parameters['Io'],trunc_val))
	print ('iip3_dbm = ',num_trunc(extracted_parameters['iip3_dbm'],trunc_val))
	print ('s11_db   = ',num_trunc(extracted_parameters['s11_db'],trunc_val))
	print ('Gain_db  = ',num_trunc(extracted_parameters['gain_db'],trunc_val))
	print ('NF_db    = ',num_trunc(extracted_parameters['nf_db'],trunc_val))
	print ('gm1      = ',num_trunc(extracted_parameters['gm1'],trunc_val))
	print ('vgs-vt   = ',num_trunc(extracted_parameters['vg']-extracted_parameters['vs']-extracted_parameters['vt'],trunc_val))
	print ('vdsat    = ',num_trunc(extracted_parameters['vdsat'],trunc_val))
	
"""
#-----------------------------------------------------------------------------------------------
# Printing the sensitivity of the parameters
def print_sensitivity(circuit_parameters_sensitivity):
	print ('\n____________________________________________________________________')
	print ('-------------------------Sensitivity --------------------------\n')
	for param_name in circuit_parameters_sensitivity:
		print('\n----------- ',param_name,' -----------\n')
		for categ in circuit_parameters_sensitivity[param_name]:
			print(categ,'\t: ',num_trunc(circuit_parameters_sensitivity[param_name][categ],trunc_val))
"""


#===========================================================================================================================