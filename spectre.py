#===========================================================================================================================
"""
Name				: Pyneni Roopesh
Roll Number			: EE18B028
File Name			: spectre.py
File Description 	: This file will contain the functions to write, run, and read from the spectre files

Functions structure in this file:
	--> valueName_to_value
	--> valueE_to_value
	--> extract_file
	--> extract_basic_parameters
		--> extract_dc_param
		--> extract_ac_param
		--> extract_sp_param
		--> extract_noise_param
	--> calculate_iip3_single_point
	--> calculate_iip3_multiple_points
		--> calculate_slope
		--> calculate_best_iip3_point
		--> check_freq
		--> extract_vout_magnitude
		--> extract_vout

	--> print_param
	--> dict_convert
	--> write_circuit_parameters
	--> write_MOS_parameters
	--> write_simulation_parameters
	--> write_tcsh_file

	--> write_extract
		--> write_extract_basic
		--> write_extract_iip3
		--> run_file


	
"""
#===========================================================================================================================
from typing import no_type_check
import numpy as np
import fileinput
import os


"""
====================================================================================================================================================================================
------------------------------------------------------------ EXTRACTION FUNCTION ---------------------------------------------------------------------------------------------------
"""

#===========================================================================================================================================================
#------------------------------------------------------ Character to Real Number Functions -----------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as a string to a floating point value 
# Input: Value of the number in string format 	
# Output: Value of the number in float
def valueName_to_value(value_name):

	# Checking if the last character of array is a string
	if value_name[-1].isalpha()==0:
		val=float(value_name)
		return val
	
	# Checking if the last character of array is a string
	if (value_name[-1]=='G' and value_name[-2]=='E') or (value_name[-1]=='g' and value_name[-2]=='e'):
		val=float(value_name[:-3])*1e6
		return val
		
	# Extracting the numerical part of the number 
	val=float(value_name[:-1])
	
	# Extracting the character that denotes the units ( i.e, millt, micro, nano, etc)
	mult_name=value_name[-1]
	mult=1.0
	
	# Calculating the value of the unit
	if mult_name=='M' or mult_name=='m':
		mult=1e-3
	elif mult_name=='U' or mult_name=='u':
		mult=1e-6
	elif mult_name=='N' or mult_name=='n':
		mult=1e-9
	elif mult_name=='P' or mult_name=='p':
		mult=1e-12
	elif mult_name=='F' or mult_name=='f':
		mult=1e-15
	elif mult_name=='G' or mult_name=='g':
		mult=1e9
	else:
		mult=1.0
		
	val=val*mult
	return val
	
#---------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as 10e1, 1.5e-2 to a floating point value 
# Input: Value of the number in string format 	
# Output: Value of the number in float
def valueE_to_value(value_name):
    
    # Extracting the number before and after e
    if 'e' in value_name:
        num1=float(value_name.split('e')[0])
        num2=float(value_name.split('e')[1])
        
        # Calculating the final number
        num=num1*(10**num2)
    
    else:
        num=float(value_name)
    
    return num




#===========================================================================================================================================================
#--------------------------------------------------------- Other File Extraction Functions -----------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
# Inputs: file name
# Output: array of lines
def extract_file(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()
	return lines


#===========================================================================================================================================================
#------------------------------------------------------ Basic File Extraction Functions --------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting the DC from the file
# Inputs: Optimization_input_parameters
# Output: Dictionary with all the parameters
def extract_dc_param(optimization_input_parameters):

	# Getting the filename
	file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/dc.out'
	lines=extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values from the required line
	extracted_parameters['vd']=valueE_to_value(lines[1])
	extracted_parameters['vg']=valueE_to_value(lines[2])
	extracted_parameters['vs']=valueE_to_value(lines[3])

	extracted_parameters['i_source']=np.absolute(valueE_to_value(lines[4]))
	extracted_parameters['v_source']=np.absolute(valueE_to_value(lines[5]))
	extracted_parameters['p_source']=extracted_parameters['i_source']*extracted_parameters['v_source']

	extracted_parameters['Io']=valueE_to_value(lines[6])
	extracted_parameters['gm1']=valueE_to_value(lines[7])
	extracted_parameters['gds1']=valueE_to_value(lines[8])
	extracted_parameters['vt']=valueE_to_value(lines[9])
	extracted_parameters['vdsat']=valueE_to_value(lines[10])

	extracted_parameters['cgd1']=np.absolute(valueE_to_value(lines[11]))
	extracted_parameters['cgs1']=np.absolute(valueE_to_value(lines[12]))

	return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting the AC from the file
# Inputs: Optimization_input_parameters
# Output: Dictionary with all the parameters
def extract_ac_param(optimization_input_parameters):

	# Getting the filename
	file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/ac.out'
	lines=extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values frim the required line
	extracted_parameters['gain_db']=valueE_to_value(lines[4])
	extracted_parameters['freq']=valueE_to_value(lines[0])

	return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting the SP from the file
# Inputs: Optimization_input_parameters
# Output: Dictionary with all the parameters
def extract_sp_param(optimization_input_parameters):

	# Getting the filename
	file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/sp.out'
	lines=extract_file(file_name)
	
	extracted_parameters={}
	
	# Skipping the first few lines
	lines=lines[12:]
	lines=lines[0].split()

	# Extracting the value from the required line
	num_char=lines[1].split(',')[0]
	extracted_parameters['s11_db']=valueE_to_value(num_char)

	return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting the Noise from the file
# Inputs: Optimization_input_parameters
# Output: Dictionary with all the parameters
def extract_noise_param(optimization_input_parameters):

	# Getting the filename
	file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/noise.out'
	lines=extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the value from the required line
	extracted_parameters['nf_db']=valueE_to_value(lines[1])

	return extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------
# Extracting all the output parameters from chi file
# Inputs: optimization_input parameters
# Outputs: output parameters dictionary 
def extract_basic_parameters(optimization_input_parameters):
	
	# Extracting the outputs 
	extracted_parameters_dc=extract_dc_param(optimization_input_parameters)
	extracted_parameters_ac=extract_ac_param(optimization_input_parameters)
	extracted_parameters_sp=extract_sp_param(optimization_input_parameters)
	extracted_parameters_noise=extract_noise_param(optimization_input_parameters)
	
	# Storing the outputs in a single dictionary
	extracted_parameters={}
	
	for param_name in extracted_parameters_dc:
		extracted_parameters[param_name]=extracted_parameters_dc[param_name]
	for param_name in extracted_parameters_ac:
		extracted_parameters[param_name]=extracted_parameters_ac[param_name]
	for param_name in extracted_parameters_sp:
		extracted_parameters[param_name]=extracted_parameters_sp[param_name]
	for param_name in extracted_parameters_noise:
		extracted_parameters[param_name]=extracted_parameters_noise[param_name]
	
	return extracted_parameters


#===========================================================================================================================================================
#------------------------------------------------------------ IIP3 Extraction Functions --------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the IIP3 from a single point
# Inputs: Vout_fund, Vout_im3, pin
# Output: IIP3
def calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin):

	# Calculating values in log scale
	vout_fund_log=20*np.log10(vout_fund_mag)
	vout_im3_log=20*np.log10(vout_im3_mag)

	# Calculating iip3
	iip3=pin+(0.5*(vout_fund_log-vout_im3_log))

	return iip3

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the IIP3 after extraction of Vout data
# Inputs: Optimization_input_parameters, Vout_fund, Vout_im3, pin
# Output: IIP3
def calculate_iip3_multiple_points(optimization_input_parameters,vout_fund_mag,vout_im3_mag,pin):

	# Calculating values in log scale
	vout_fund_log=20*np.log10(vout_fund_mag)
	vout_im3_log=20*np.log10(vout_im3_mag)

	n_pin=optimization_input_parameters['simulation']['pin_points']
	n_points=optimization_input_parameters['simulation']['iip3_calc_points']
	
	# Creating arrays for slopes and y-intercepts of fundamental and im3 components
	fund_slope=np.zeros(n_pin+1-n_points,dtype=float)
	fund_intercept=np.zeros(n_pin+1-n_points,dtype=float)
	im3_slope=np.zeros(n_pin+1-n_points,dtype=float)
	im3_intercept=np.zeros(n_pin+1-n_points,dtype=float)

	# Calculating the slopes and y-intercepts
	for i in range(n_pin+1-n_points):
		fund_slope[i],fund_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_fund_log[i:i+n_points-1])
		im3_slope[i],im3_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_im3_log[i:i+n_points-1])
	
	# Finding the best points for iip3 calculation
	best_point=calculate_best_iip3_point(fund_slope,im3_slope)
	
	# Calculating the iip3 given the slope and y-intercept of fundamental and im3
	iip3=(im3_intercept[best_point]-fund_intercept[best_point])/(fund_slope[best_point]-im3_slope[best_point])

	return iip3

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the slope and y-intercept
# Inputs: x and y coordinates of the points
# Output: slope, y-intercept
def calculate_slope(x,y):
	A = np.vstack([x, np.ones(len(x))]).T
	m, c = np.linalg.lstsq(A, y, rcond=None)[0]
	return m,c

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the point with slope closest to 1dB/dB for fund and 3dB/dB for im3
# Inputs: Slope of fundamental and im3
# Output: Location of the best point
def calculate_best_iip3_point(fund_slope,im3_slope):
	
	# Getting the length of the array
	l=len(fund_slope)

	# Choosing the best point as the first point
	best_point=0
	best_error=(fund_slope[0]-1)**2+(im3_slope[0]-3)**2

	for i in range(1,l):
		error=(fund_slope[i]-1)**2+(im3_slope[i]-3)**2
		if error<best_error:	# Checking if the current point is better than the best point
			best_point=i
			best_error=error
	return best_point

#---------------------------------------------------------------------------------------------------------------------------	
# Checks if the frequency is within range ( within (target-error,target+error) )
# Inputs: Test Frequency, Target Frequency, Error
# Output: 1 if Yes and 0 if No
def check_freq(f_test,f_target,f_error):
	if f_test<f_target+f_error and f_test>f_target-f_error:
		return 1
	else:
		return 0

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting Vout magnitude of fundamental and im3 from file ( for hb_sweep )
# Inputs: Filename, Optimization Input Parameters
# Output: Magnitude of Vout at fundamental and im3
def extract_vout_magnitude(file_name,optimization_input_parameters):

	lines=extract_file(file_name)

	fund_1=optimization_input_parameters['simulation']['parameters_list']['fund_1']
	fund_2=optimization_input_parameters['simulation']['parameters_list']['fund_2']

	f_im3=2*fund_2-fund_1
	f_error=(fund_2-fund_1)/100

	flag=0
	flag_fun=0
	flag_im3=0
	flag_test=0

	while 1:
		if len(lines[0].split())<2:
			lines=lines[1:]
		
		elif 'freq' in lines[0].split()[0] and flag==0:
			flag=1
			lines=lines[1:]
		
		elif 'freq' in lines[0].split()[0] and flag==1:
			if flag_fun==0 and check_freq(float(lines[0].split()[1]),fund_2,f_error)==1 :
				
				#Extracting Vout for fundamental
				flag_test=1
				while flag_test==1:
					if 'Vout' in lines[0].split()[0]:
						flag_test=0
						vout_fund=extract_vout(lines[0])
					else:
						lines=lines[1:]
				flag_fun=1
			
			elif flag_im3==0 and check_freq(float(lines[0].split()[1]),f_im3,f_error)==1 :
				
				#Extracting Vout for im3
				flag_test=1
				while flag_test==1:
					if 'Vout' in lines[0].split()[0]:
						flag_test=0
						vout_im3=extract_vout(lines[0])
					else:
						lines=lines[1:]
				flag_im3=1
			lines=lines[1:]
			
			if flag_fun==1 and flag_im3==1:
				break
		
		else:
			lines=lines[1:]

	return vout_fund,vout_im3

#---------------------------------------------------------------------------------------------------------------------------	
# Extracts Vout_magnitude from hb,pss file line
# Inputs: Line
# Output: Vout_Magnitude
def extract_vout(lines):
	
	# Extracting Vout Magnitude
	lines=lines.split()
	char_r=lines[1].split('(')[1]
	char_i=lines[2].split(')')[0]

	# Converting string to floating point value
	vout_r=valueE_to_value(char_r)
	vout_i=valueE_to_value(char_i)
	
	# Calculating the magnitude of the output
	vout_mag=np.sqrt(vout_r*vout_r+vout_i*vout_i)

	return vout_mag


#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ FILE WRITE FUNCTIONS --------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Command that returns the string that has to be printed in the .scs file
# Inputs  : Name of the parameter, Value of the parameter
# Outputs : String to be printed
def print_param(param_var,val):
	return "parameters "+param_var+'='+str(val)+'\n'

#-----------------------------------------------------------------      
# Function that converts input parameter dictionary to writing dictionary
# Inputs  : Circuit Parameters Dictionary, Optimization Input Parameters
# Outputs : The dictionary containing the parameters to be written to the .scs file
def dict_convert(circuit_parameters,optimization_input_parameters):
	write_dict={}
	
	# param_names in write_dict will contain the name of the parameters as it is written in the .scs file
	for param_name in optimization_input_parameters['simulation']['cir_writing_dict']:
		write_dict[param_name]=circuit_parameters[optimization_input_parameters['simulation']['cir_writing_dict'][param_name]]

	# Checking if we have TSMC Resistors
	if optimization_input_parameters['simulation']['basic_circuit']=='basic_parameters_tsmc_65':
		write_dict['Resb_L'],write_dict['Resb_W']=get_TSMC_resistor(circuit_parameters['Rb'])
		write_dict['Resd_L'],write_dict['Resd_W']=get_TSMC_resistor(circuit_parameters['Rd'])
		write_dict['Resbias_L'],write_dict['Resbias_W']=get_TSMC_resistor(circuit_parameters['Rbias'])
	
	# Calculating the number of fingers
	n_finger=int(circuit_parameters['W']/optimization_input_parameters['simulation']['w_finger_max'])+1
	write_dict['n_finger']=n_finger
	
	return write_dict

#-----------------------------------------------------------------      
# Function that converts resistance to length and width
# Inputs  : resistance
# Outputs : length, width
def get_TSMC_resistor(resistance):
	sheet_resistance=124.45
	#L_min=0.8e-6
	W_min=0.4e-6
	dW=0.0691e-6

	"""
	if resistance<sheet_resistance:
		length=L_min
		width=L_min*sheet_resistance/resistance
	else:
		width=W_min
		length=W_min*resistance/sheet_resistance
	"""
	width=W_min-dW
	length=width*resistance/sheet_resistance
	
	return length,width
            
#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_circuit_parameters(circuit_parameters,optimization_input_parameters):
	
	# We will convert the circuit parameters to write_dict
	write_dict=dict_convert(circuit_parameters,optimization_input_parameters)
	
	# Getting the filenames
	filename1=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/circ.scs'
	filename2=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'/circ.scs'

	# We will write the new values to the Basic Circuit
	f=open(filename1,'r+')
	s=''
	for line in fileinput.input(filename1):
		for param_name in write_dict:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,print_param(param_name,write_dict[param_name]))	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

	# We will write the new values to the IIP3 Circuit
	f=open(filename2,'r+')
	s=''
	for line in fileinput.input(filename2):
		for param_name in write_dict:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,print_param(param_name,write_dict[param_name]))	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that adds MOSFET Parameters to the netlist
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def write_MOS_parameters(optimization_input_parameters):
	
	# write_dict will contain the mosfet values
	write_dict={
		'len':optimization_input_parameters['MOS']['Lmin'],
		'v_dd':optimization_input_parameters['MOS']['Vdd'],
	}
	process_corner=optimization_input_parameters['simulation']['process_corner']

	# Getting the filenames
	filename1=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/circ.scs'
	filename2=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'/circ.scs'

	# Writing the MOS Parameters to Basic File
	f=open(filename1,'r+')
	s=''
	write_check=1
	include_check=0

	# Replacing the lines of .scs file
	for line in fileinput.input(filename1):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+optimization_input_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,print_param(param_name,write_dict[param_name]))

		if write_check==1:
			s=s+line

	f.truncate(0)
	f.write(s)
	f.close()

	# Writing the MOS Parameters to IIP3 File
	f=open(filename2,'r+')
	s=''
	write_check=1
	include_check=0

	# Replacing the lines of .scs file
	for line in fileinput.input(filename2):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+optimization_input_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,print_param(param_name,write_dict[param_name]))
		
		if 'hb_test' in line and 'errpreset=conservative' in line and optimization_input_parameters['simulation']['conservative']=='NO':
			line_split=line.split()
			line=''
			for word in line_split:
				if 'errpreset=conservative' not in word:
					line=line+word+' '
			line=line+'\n'
		
		elif 'hb_test' in line and 'errpreset=conservative' not in line and optimization_input_parameters['simulation']['conservative']=='YES':
			line=line[:-1]+' errpreset=conservative \n'

		if write_check==1:
			s=s+line
	

	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that adds Simulation Parameters
# Inputs  : Optimization Input Parameters
# Outputs : NONE
def write_simulation_parameters(optimization_input_parameters):
	
	# Adding simulation_parameters to write_dict
	write_dict={}
	for param_name in optimization_input_parameters['simulation']['parameters_list']:
		write_dict[param_name]=optimization_input_parameters['simulation']['parameters_list'][param_name]
	process_corner=optimization_input_parameters['simulation']['process_corner']

	# Getting the filenames
	filename1=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'/circ.scs'
	filename2=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'/circ.scs'

	# Writing the simulation parameters to Basic File
	f=open(filename1,'r+')
	s=''
	write_check=1
	include_check=0
	
	# Replacing the lines of .scs file
	for line in fileinput.input(filename1):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+optimization_input_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,print_param(param_name,write_dict[param_name]))
		
		if write_check==1:
			s=s+line

	f.truncate(0)
	f.write(s)
	f.close()

	# Writing the simulation parameters to IIP3 File
	f=open(filename2,'r+')
	s=''
	write_check=1
	include_check=0

	# Replacing the lines of .scs file
	for line in fileinput.input(filename2):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+optimization_input_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,print_param(param_name,write_dict[param_name]))
				
		if 'hb_test' in line and 'errpreset=conservative' in line and optimization_input_parameters['simulation']['conservative']=='NO':
			line_split=line.split()
			line=''
			for word in line_split:
				if 'errpreset=conservative' not in word:
					line=line+word+' '
			line=line+'\n'
		
		elif 'hb_test' in line and 'errpreset=conservative' not in line and optimization_input_parameters['simulation']['conservative']=='YES':
			line=line[:-1]+' errpreset=conservative \n'
			
		
		if write_check==1:
			s=s+line

	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that modifies tcsh file
# Inputs  : Optimization_Input_Parameters
# Outputs : NONE
def write_tcsh_file(optimization_input_parameters,optimiztion_type):
	
	filename=optimization_input_parameters['simulation']['tcsh']
	f=open(filename,'r+')
	s=''
	
	s='#tcsh\n'
	s=s+'source ~/.cshrc\n'
	
	if optimiztion_type=='basic':
		s=s+'cd '+optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['basic_circuit']+'\n'
	else:
		s=s+'cd '+optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'\n'
	
	s=s+'spectre circ.scs =log circ_log.txt\n'
	s=s+'exit'
	
	f.truncate(0)
	f.write(s)
	f.close()

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ SPECTRE RUNNING FUNCTIONS ---------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------	
# This function will run the shell commands to run Spectre
# Inputs  : Optimization Input Parameters
# Outputs : NONE

def run_file(optimization_input_parameters):
	os.system('cd /home/ee18b028/cadence_project')
	os.system('tcsh '+optimization_input_parameters['simulation']['tcsh'])	# This is the command to run the spectre file
	
#-----------------------------------------------------------------------------------------------
# This function will perform simulation for Basic Parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract_basic(optimization_input_parameters):
	
	# Writing the simulation parameters
	write_simulation_parameters(optimization_input_parameters)

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(optimization_input_parameters,'basic')

	# Running netlist file
	run_file(optimization_input_parameters)

	# Extracting the Basic Parameters
	basic_extracted_parameters=extract_basic_parameters(optimization_input_parameters)
	
	return basic_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will perform simulation for IIP3 Parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract_iip3(optimization_input_parameters):
	
	if optimization_input_parameters['simulation']['iip3_type']=='basic':
		
		pin=optimization_input_parameters['simulation']['pin_fixed']
		optimization_input_parameters['simulation']['parameters_list']['pin']=pin
		
		# Writing the simulation parameters
		write_simulation_parameters(optimization_input_parameters)

		# Writing the tcsh file for Basic Analysis
		write_tcsh_file(optimization_input_parameters,'iip3')

		# Running netlist file
		run_file(optimization_input_parameters)

		# Extracting Vout Magnitude
		file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'/circ.raw/hb_test.fd.qpss_hb'
		vout_fund_mag,vout_im3_mag=extract_vout_magnitude(file_name,optimization_input_parameters)

		# Calculating the iip3
		iip3=calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin)

	else:
		pin_start=optimization_input_parameters['simulation']['pin_start']
		pin_stop=optimization_input_parameters['simulation']['pin_stop']
		pin_points=optimization_input_parameters['simulation']['pin_points']

		pin=np.linspace(pin_start,pin_stop,pin_points)
		
		vout_fund_mag=np.zeros(pin_points,dtype=float)
		vout_im3_mag=np.zeros(pin_points,dtype=float)

		for i in range(pin_points):
		
			optimization_input_parameters['simulation']['parameters_list']['pin']=pin[i]
			
			# Writing the simulation parameters
			write_simulation_parameters(optimization_input_parameters)

			# Writing the tcsh file for Basic Analysis
			write_tcsh_file(optimization_input_parameters,'iip3')

			# Running netlist file
			run_file(optimization_input_parameters)

			# Extracting Vout Magnitude
			file_name=optimization_input_parameters['simulation']['directory']+optimization_input_parameters['simulation']['iip3_circuit']+'/circ.raw/hb_test.fd.qpss_hb'
			vout_fund_mag[i],vout_im3_mag[i]=extract_vout_magnitude(file_name,optimization_input_parameters)

		iip3=calculate_iip3_multiple_points(optimization_input_parameters,vout_fund_mag,vout_im3_mag,pin)

	iip3_extracted_parameters={'iip3_dbm':iip3}
	
	return iip3_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract(circuit_parameters,optimization_input_parameters):
	
	# Writing to netlist file
	write_circuit_parameters(circuit_parameters,optimization_input_parameters)

	# Extracting the Basic Parameters
	basic_extracted_parameters=write_extract_basic(optimization_input_parameters)

	# Extracting the IIP3 Parameters
	iip3_extracted_parameters=write_extract_iip3(optimization_input_parameters)

	# Extracting Parameters from output files
	extracted_parameters=basic_extracted_parameters.copy()
	for param_name in iip3_extracted_parameters:
		extracted_parameters[param_name]=iip3_extracted_parameters[param_name]
	
	return extracted_parameters

#===========================================================================================================================
