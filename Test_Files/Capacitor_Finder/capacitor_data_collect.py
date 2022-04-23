#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file is used to store the value of Q and L of TSMC Inductor by sweeping various inductor parameters
"""

#===========================================================================================================================
import numpy as np
import fileinput
import os
from matplotlib import pylab
from pylab import *


"""
============================================================================================================================
------------------------------------------ EXTRACTION FUNCTION -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as 10e1, 1.5e-2 to a floating point value 
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

#---------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
def extract_file(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()
	return lines

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
# Extracting the HB from the file
# Inputs: Optimization_input_parameters
# Output: resistance_dict,distortion,symmetry
def extract_hb_param(filename,freq,i_cur):

	# Getting the lines in the filename
	lines=extract_file(filename)
	capacitance=0
	
	# Extracting the data from the HB Raw File
	flag=0
	for line in lines:
		
		# If flag=1, it means that the "freq" line is already over and we need to extract the vout values.
		if flag==1:
			
			# Extracting R1a node values
			if 'Rn1' in line:
				char_real=(line.split()[1])[1:]
				char_img=(line.split()[2])[:-1]
				vout_a_real_array=-1*float(char_img)
				vout_a_img_array=float(char_real)
			
			# Extracting R1b node values
			elif 'Rn2' in line:
				char_real=(line.split()[1])[1:]
				char_img=(line.split()[2])[:-1]
				vout_b_real_array=-1*float(char_img)
				vout_b_img_array=float(char_real)
				break
			else:
				continue
		
		if 'freq' in line:
			if 'sweep' in line:
				continue
			frequency=float(line.split()[1])
			if check_freq(frequency,freq,freq*0.001)==1:
				flag=1

	voltage_diff=np.sqrt((vout_a_real_array-vout_b_real_array)**2+(vout_a_img_array-vout_b_img_array)**2)
	capacitance=i_cur/(2*np.pi*freq*voltage_diff)
	
	return capacitance

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the AC from the file
def extract_ac_param(filename):

	# Getting the filename
	lines=extract_file(filename)

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values frim the required line
	capacitance=valueE_to_value(lines[6])
	
	return capacitance

"""
============================================================================================================================
------------------------------------------ FILE WRITE FUNCTIONS ------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_circuit_parameters(filename,circuit_parameters):
	
	# We will write the new values to the netlist file
	f=open(filename,'r+')
	s=''
	for line in fileinput.input(filename):
		for param in circuit_parameters:
			if 'parameters '+str(param) in line:										# Checking for a particular parameter in the .scs file
				line='parameters '+str(param)+'='+str(circuit_parameters[param])+'\n'	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that modifies tcsh file
# Inputs  : Optimization_Input_Parameters
# Outputs : NONE
def write_tcsh_file(file_directory):
	
	filename='/home/ee18b028/cadence_project/Spectre_Run/spectre_run.tcsh'
	f=open(filename,'r+')
	s=''
	
	s='#tcsh\n'
	s=s+'source ~/.cshrc\n'
	s=s+'cd '+file_directory+'\n'
	s=s+'spectre circ.scs =log circ_log.txt\n'
	#s=s+'spectre circ.scs\n'
	s=s+'exit'
	
	f.truncate(0)
	f.write(s)
	f.close()


"""
============================================================================================================================
------------------------------------------ SPECTRE RUNNING FUNCTIONS -------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------	
# This function will run the shell commands to run Spectre
def run_file():
	os.system('tcsh /home/ee18b028/cadence_project/Spectre_Run/spectre_run.tcsh')	# This is the command to run the spectre file
	
#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters
def write_extract(file_directory,circuit_parameters):

	# Getting the filenames
	filename_w=file_directory+'/circ.scs'
	filename_r=file_directory+'/circ.raw/hb_test.fd.pss_hb'
	#filename_r=file_directory+'/ac.out'

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(file_directory)

	# Writing to netlist file
	write_circuit_parameters(filename_w,circuit_parameters)

	# Running netlist file
	run_file()

	# Extracting the HB Parameters
	freq=circuit_parameters['fund_1']
	i_cur=circuit_parameters['i_sin']
	
	capacitance=extract_hb_param(filename_r,freq,i_cur)
	#capacitance=extract_ac_param(filename_r)
	
	return capacitance


"""
============================================================================================================================
------------------------------------------ SWEEPING CIRCUIT PARAMETERS -----------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the Q and L of the inductor for various circuit parameters
def Parameter_Sweep(file_directory_netlist,file_directory_output):

	circuit_parameters={
		'len':0.2e-6,
		'wid':0.2e-6,
		'i_sin':1e-6,
		'm_factor':1.0,
		'fund_1':1e9
	}

	len_array=np.linspace(5e-6,100e-6,20)
	wid_array=np.linspace(5e-6,10e-6,2)
	
	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_output):
		os.makedirs(file_directory_output)

	# Opening the file
	filename_csv=file_directory_output+'capacitor_sweep.csv'
	f=open(filename_csv,'w')
	f.write('Width,Length,Capacitance\n')

	for wid in wid_array:
		for len in len_array:
			
						
			# Storing the circuit parameters
			circuit_parameters={
				'len':len,
				'wid':wid,
				'i_sin':1e-6,
				'm_factor':1.0,
				'fund_1':1e9
			}
			
			print(circuit_parameters)
			
			# Running spectre
			capacitance=write_extract(file_directory_netlist,circuit_parameters)

			# Writing the results
			f.write(str(wid)+',')
			f.write(str(len)+',')
			f.write(str(capacitance)+'\n')

	f.close()
		
	
"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/Test_Circuits/Capacitor/capacitor_test_1'

start=datetime.datetime.now()
print(start)

# Code to do capacitor sweep data storage
file_directory_output='/home/ee18b028/Optimization/Simulation_Results/Capacitor/Sweep2/'
Parameter_Sweep(file_directory,file_directory_output)

end=datetime.datetime.now()
print(start)
print(end)
print(end-start)

