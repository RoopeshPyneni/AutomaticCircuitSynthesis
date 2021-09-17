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
	resistance=valueE_to_value(lines[2])

	return resistance

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ FILE WRITE FUNCTIONS --------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_resistor_name(filename,resistor_name):
	
	# We will write the new values to the Basic Circuit
	f=open(filename,'r+')
	s=''
	for line in fileinput.input(filename):
		if "R1" in line:	# Checking for a particular parameter in the .scs file
			line=line.split()
			if line[0]=='R1':
				line[3]=resistor_name
			newline=''
			for word in line:
				newline=newline+word+' '
			line=newline+'\n'
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_circuit_parameters(filename,len,wid,temp):
	
	# Creating a dictionary of the parameters
	write_dict={'len':len,'wid':wid,'cir_temp':temp}

	# We will write the new values to the netlist file
	f=open(filename,'r+')
	s=''
	for line in fileinput.input(filename):
		for param in write_dict:
			if 'parameters '+str(param) in line:								# Checking for a particular parameter in the .scs file
				line='parameters '+str(param)+'='+str(write_dict[param])+'\n'	# Replacing the parameter in the .scs file
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
	s=s+'cd /home/ee18b028/cadence_project/test/resistor_test\n'
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
	os.system('tcsh /home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh')	# This is the command to run the spectre file
	
#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract(filename_w,filename_e,len,wid,temp):
	
	# Writing the tcsh file for Basic Analysis
	write_tcsh_file()

	# Writing to netlist file
	write_circuit_parameters(filename_w,len,wid,temp)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	resistance=extract_dc_param(filename_e)
	
	return resistance

#===========================================================================================================================

"""
====================================================================================================================================================================================
------------------------------------------------------------ PLOTS -----------------------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will plot all the graphs
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def plot_resistance(file_directory_plot,resistance_array,temp_array,len_array,wid_array):
	
	# First, we will plot R vs T for different l 
	print('Start Plot vs Temperature')
	file_directory=file_directory_plot+'/Temperature'
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)

	for k in range(len(wid_array)):
		figure()
		for j in range(len(len_array)):
			plot(temp_array,resistance_array[:,j,k],label='Length = '+str(len_array[j]))
		xlabel('Temperature')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/wid_'+str(wid_array[k])+'.pdf')
		close()

	# Second, we will plot R vs l for different w
	print('Start Plot vs Length')
	file_directory=file_directory_plot+'/Length'
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)

	for i in range(len(temp_array)):
		figure()
		for k in range(len(wid_array)):
			semilogx(len_array,resistance_array[i,:,k],label='Width = '+str(wid_array[k]))
		xlabel('Length')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/temp_'+str(temp_array[i])+'.pdf')
		close()

	# Third, we will plot R vs w for different l
	print('Start Plot vs Width')
	file_directory=file_directory_plot+'/Width'
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)

	for i in range(len(temp_array)):
		figure()
		for j in range(len(len_array)):
			semilogx(wid_array,resistance_array[i,j,:],label='Length = '+str(len_array[j]))
		xlabel('Width')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/temp_'+str(temp_array[i])+'.pdf')
		close()

	# Finally, we will plot R vs l/w for different w
	print('Start Plot vs L by W')
	file_directory=file_directory_plot+'/L_by_W'
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)

	for i in range(len(temp_array)):
		figure()
		for k in range(len(wid_array)):
			semilogx(len_array/wid_array[k],resistance_array[i,:,k],label='Width = '+str(wid_array[k]))
		xlabel('Length/Width')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/temp_'+str(temp_array[i])+'.pdf')
		close()



"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/test/resistor_test'
filename_w=file_directory+'/circ.scs'
filename_e=file_directory+'/dc.out'

# Creating the temperature, length, and width arrays
resistor_list=['rppolywo','rppolyl','rpodwo','rpodl','rnwsti','rnwod','rnpolywo','rnpolyl','rnodwo','rnodl']
temp_array=np.linspace(-40,120,17)
lw_start=60e-9
lw_end=60e-6
len_array=lw_start*np.logspace(0,np.log10(lw_end/lw_start),10)
wid_array=len_array

# Running the code 
for resistor in resistor_list:
	print('Resistor name is : ',resistor)
	write_resistor_name(filename_w,resistor)

	file_directory_plot='/home/ee18b028/Optimization/Simulation_Results/Resistance/'+resistor
	if not os.path.exists(file_directory_plot):
		os.makedirs(file_directory_plot)

	l_temp=len(temp_array)
	l_len=len(len_array)
	l_wid=len(wid_array)

	resistance_array=np.zeros((l_temp,l_len,l_wid),dtype=float)

	for i in range(l_temp):
		for j in range(l_len):
			for k in range(l_wid):
				resistance_array[i,j,k]=write_extract(filename_w,filename_e,len_array[j],wid_array[k],temp_array[i])

	plot_resistance(file_directory_plot,resistance_array,temp_array,len_array,wid_array)


"""
write_resistor_name(filename_w,'rnpolyl')

len=1200e-9
wid=6000e-9
temp=27
resistance=write_extract(filename_w,filename_e,len,wid,temp)

print('Resistance is:',resistance)
"""


"""
List of Resistor Commands in the netlist file
R0 (net017 net016 ) rppolywo l=10u w=2u m=1 mf=(1) mismatchflag=0
R1 (net017 net016 ) rppolyl l=10u w=2u m=1 mf=(1) mismatchflag=0
R2 (net017 net016 ) rpodwo l=10u w=2u m=1 mf=(1) mismatchflag=0
R3 (net016 net028 ) rpodl l=10u w=2u m=1 mf=(1) mismatchflag=0
R4 (net016 net028 ) rnwsti l=10u w=2u mf=(1)
R5 (net016 net028 ) rnwod l=10u w=2u mf=(1)
R6 (net016 net028 ) rnpolywo l=10u w=2u m=1 mf=(1) mismatchflag=0
R7 (net016 net028 ) rnpolyl l=10u w=2u m=1 mf=(1) mismatchflag=0
R8 (net016 net028 ) rnodwo l=10u w=2u m=1 mf=(1) mismatchflag=0
R9 (net033 net028 ) rnodl l=10u w=2u m=1 mf=(1) mismatchflag=0
"""