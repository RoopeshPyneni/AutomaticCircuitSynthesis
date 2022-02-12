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

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the AC from the file
def extract_ac_param(filename):

	# Getting the filename
	lines=extract_file(filename)

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values frim the required line
	Q=valueE_to_value(lines[5])
	L=valueE_to_value(lines[6])
	
	return Q,L


"""
============================================================================================================================
------------------------------------------ FILE WRITE FUNCTIONS ------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that modifies the .scs file
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
def write_tcsh_file(file_directory):
	filename='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run_extra.tcsh'
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
	os.system('tcsh /home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run_extra.tcsh')	# This is the command to run the spectre file
	
#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters
def write_extract(file_directory,circuit_parameters):

	# Getting the filenames
	filename_w=file_directory+'/circ.scs'
	filename_r=file_directory+'/ac.out'

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(file_directory)

	# Writing to netlist file
	write_circuit_parameters(filename_w,circuit_parameters)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	Q,L=extract_ac_param(filename_r)

	return Q,L


"""
============================================================================================================================
------------------------------------------ SINGLE POINT Q AND L ------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the Q and L of the inductor for the given circuit parameters
def Inductor_Single_Point(file_directory_netlist,circuit_parameters):

	# Running spectre
	Q,L=write_extract(file_directory_netlist,circuit_parameters)

	# Storing the values in arrays
	print('Q = ',Q)
	print('L = ',L)

"""
============================================================================================================================
------------------------------------------ SWEEPING CIRCUIT PARAMETERS -----------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the Q and L of the inductor for various circuit parameters
def Inductor_Circuit_Parameter_Sweep(file_directory_netlist,file_directory_output):

	circuit_parameters={
		'L_w':10e-6,
		'L_rad':30e-6,
		'L_turn':1.25,
		'L_gdis':20e-6,
		'L_spc':3e-6
	}

	w_array=np.linspace(3e-6,30e-6,10)
	rad_array=np.linspace(15e-6,90e-6,16)
	turn_array=np.linspace(0.5,5.5,21)
	gdis_array=np.linspace(10e-6,50e-6,9)
	spc_array=np.linspace(2e-6,4e-6,3)

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_output):
		os.makedirs(file_directory_output)

	# Opening the file
	filename_csv=file_directory_output+'inductor_sweep.csv'
	f=open(filename_csv,'w')
	f.write('Width,Radius,N_Turns,gdis,spc,Q,L\n')

	for w in w_array:
		for rad in rad_array:
			for turn in turn_array:
				for gdis in gdis_array:
					for spc in spc_array:
						
						# Storing the circuit parameters
						circuit_parameters={
							'L_w':w,
							'L_rad':rad,
							'L_turn':turn,
							'L_gdis':gdis,
							'L_spc':spc
						}

						# Running spectre
						Q,L=write_extract(file_directory_netlist,circuit_parameters)

						# Writing the results
						f.write(str(w)+',')
						f.write(str(rad)+',')
						f.write(str(turn)+',')
						f.write(str(gdis)+',')
						f.write(str(spc)+',')
						f.write(str(Q)+',')
						f.write(str(L)+'\n')

	f.close()
		
	
"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/Test_Circuits/Inductor/Inductor_Test_1'

# Code to find the single point Q and L
circuit_parameters={
	'L_w':10e-6,
	'L_rad':30e-6,
	'L_turn':1.25,
	'L_gdis':20e-6,
	'L_spc':3e-6
}
Inductor_Single_Point(file_directory,circuit_parameters)


# Code to do distortion analysis for multiple
file_directory_output='/home/ee18b028/Optimization/Simulation_Results/Inductor/Sweep/'
Inductor_Circuit_Parameter_Sweep(file_directory,file_directory_output)

