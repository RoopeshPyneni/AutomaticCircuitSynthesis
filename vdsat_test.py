#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Vdsat File:
"""
#===========================================================================================================================
from typing import no_type_check
import numpy as np
import fileinput
import os
from matplotlib import pylab
from pylab import *

"""
====================================================================================================================================================================================
------------------------------------------------------------ EXTRACTION FUNCTION ---------------------------------------------------------------------------------------------------
"""

#===========================================================================================================================================================
#------------------------------------------------------ Character to Real Number Functions -----------------------------------------------------------------

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

def extract_dc_param(filename):

	lines=extract_file(filename)
	
	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values
	id=valueE_to_value(lines[4])
	vdsat=valueE_to_value(lines[6])

	return id,vdsat

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ FILE WRITE FUNCTIONS --------------------------------------------------------------------------------------------------
"""
      
#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_circuit_parameters(filename,v_g):
	
	# We will write the new values to the Basic Circuit
	f=open(filename,'r+')
	s=''
	for line in fileinput.input(filename):
		if "parameters v_g=" in line:	# Checking for a particular parameter in the .scs file
			line='parameters v_g='+str(v_g)	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that modifies tcsh file
# Inputs  : Optimization_Input_Parameters
# Outputs : NONE
def write_tcsh_file():
	
	filename='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	f=open(filename,'r+')
	s=''
	
	s='#tcsh\n'
	s=s+'source ~/.cshrc\n'
	s=s+'cd /home/ee18b028/cadence_project/lna1/vdsat_test'
	s=s+'spectre circ.scs \n'
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

def run_file():
	os.system('cd /home/ee18b028/cadence_project')
	os.system('tcsh /home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh')	# This is the command to run the spectre file


#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract(filename_w,v_g,filename_e):
	
	# Writing the tcsh file for Basic Analysis
	write_tcsh_file()

	# Writing to netlist file
	write_circuit_parameters(filename_w,v_g)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	id,vdsat=extract_dc_param(filename_e)
	
	return id,vdsat

#===========================================================================================================================

"""
====================================================================================================================================================================================
------------------------------------------------------------ PLOT FUNCTIONS --------------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def plot_id_vg(file_directory_plot,id,vg):
	
	figure()
	plot(vg,id,label='id')
	ylabel('id')
	xlabel('vg')
	title('id vs vg')
	grid()
	legend()
	savefig(file_directory_plot+'id_vs_vg.pdf')
	close()



"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

id_array=[]
vdsat_array=[]
vg_array=[]

file_directory='/home/ee18b028/cadence_project/lna1/vdsat_test'
filename_w=file_directory+'/circ.scs'
filename_e=file_directory+'/dc.out'

for vg in np.linspace(0,1,101):
	id,vdsat=write_extract(filename_w,vg,filename_e)
	id_array.append(id)
	vdsat_array.append(vdsat)
	vg_array.append(vg)

id_array=np.array(id_array)
vdsat_array=np.array(vdsat_array)
vg_array=np.array(vg_array)

file_directory_plot='/home/ee18b028/Optimization/Simulation_Results/Vdsat/'
plot_id_vg(file_directory_plot,id,vg)