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
	i_d=valueE_to_value(lines[4])
	vdsat=valueE_to_value(lines[6])

	return i_d,vdsat

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
			line='parameters v_g='+str(v_g)+'\n'	# Replacing the parameter in the .scs file
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
	s=s+'cd /home/ee18b028/cadence_project/lna1/vdsat_test\n'
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
	#os.system('cd /home/ee18b028/cadence_project')
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
	i_d,vdsat=extract_dc_param(filename_e)
	
	return i_d,vdsat

#===========================================================================================================================

"""
====================================================================================================================================================================================
------------------------------------------------------------ DIFFERENTIATION -------------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# 
# Inputs  : 
# Outputs : 

def differentiation_array(arr_x,arr_y):
	
	arr_y1=np.zeros(len(arr_y)-1,dtype=float)
	arr_y2=np.zeros(len(arr_y)-2,dtype=float)
	arr_y3=np.zeros(len(arr_y)-3,dtype=float)

	for i in range(len(arr_y1)):
		arr_y1[i]=(arr_y[i+1]-arr_y[i])/(arr_x[i+1]-arr_x[i])
	
	for i in range(len(arr_y2)):
		arr_y2[i]=(arr_y1[i+1]-arr_y1[i])/(arr_x[i+1]-arr_x[i])
	
	for i in range(len(arr_y3)):
		arr_y3[i]=(arr_y2[i+1]-arr_y2[i])/(arr_x[i+1]-arr_x[i])

	return arr_x[:-3],arr_y[:-3],arr_y1[:-2],arr_y2[:-1],arr_y3

"""
====================================================================================================================================================================================
------------------------------------------------------------ PLOT FUNCTIONS --------------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# 
# Inputs  : 
# Outputs : 

def plot_curves(file_directory_plot,vg,i_d,gm1,gm2,gm3):
	
	# ID vs VG
	figure()
	plot(vg,i_d)
	ylabel('id')
	xlabel('vg')
	title('id vs vg')
	grid()
	savefig(file_directory_plot+'/id_vs_vg.pdf')
	close()

	# GM1 vs VG
	figure()
	plot(vg,gm1)
	ylabel('gm1')
	xlabel('vg')
	title('gm1 vs vg')
	grid()
	savefig(file_directory_plot+'/gm1_vs_vg.pdf')
	close()

	# GM2 vs VG
	figure()
	plot(vg,gm2)
	ylabel('gm2')
	xlabel('vg')
	title('gm2 vs vg')
	grid()
	savefig(file_directory_plot+'/gm2_vs_vg.pdf')
	close()

	# GM3 vs VG
	figure()
	plot(vg,gm3)
	ylabel('gm3')
	xlabel('vg')
	title('gm3 vs vg')
	grid()
	savefig(file_directory_plot+'/gm3_vs_vg.pdf')
	close()

	# All vs VG
	figure()
	plot(vg,i_d,label='id')
	plot(vg,gm1,label='gm1')
	plot(vg,gm2,label='gm2')
	plot(vg,gm3,label='gm3')
	ylabel('All Parameters')
	xlabel('vg')
	title('All Parameters vs vg')
	grid()
	legend()
	savefig(file_directory_plot+'/all_vs_vg.pdf')
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
	i_d,vdsat=write_extract(filename_w,vg,filename_e)
	id_array.append(i_d)
	vdsat_array.append(vdsat)
	vg_array.append(vg)

id_array=np.array(id_array)
vdsat_array=np.array(vdsat_array)
vg_array=np.array(vg_array)

file_directory_plot='/home/ee18b028/Optimization/Simulation_Results/Vdsat'
if not os.path.exists(file_directory_plot):
	os.makedirs(file_directory_plot)

vg_array,id_array,gm1_array,gm2_array,gm3_array=differentiation_array(vg_array,id_array)
vdsat_array=vdsat_array[:-3]

plot_curves(file_directory_plot,vg_array,id_array,gm1_array,gm2_array,gm3_array)

vdsat_best=0
i_best=0
vg_best=0
id_best=9
for i in range(len(vg_array)-1):
	if gm3_array[i]>0 and gm3_array[i+1]<=0:
		vdsat_best=vdsat_array[i]
		i_best=i
		vg_best=vg_array[i]
		id_best=id_array[i]
		break

print('Vdsat Best : ',vdsat_best)
print('i     Best : ',i_best)
print('Vg    Best : ',vg_best)
print('Id    Best : ',id_best)
