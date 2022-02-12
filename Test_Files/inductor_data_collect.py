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
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE DISTORTION -------------------------------------------------------------------------------------------
"""

""""
#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the distortion for all the resistors over extreme sizes
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
# UPTO DATE
def MOS_Resistor_Distortion(file_directory_netlist,resistor_dict,file_directory_output):

	# Creating current array
	current_array=np.logspace(-6,-3,7)
	vout_fund_array=np.zeros(len(current_array),dtype=float)
	vout_harm_array=np.zeros(len(current_array),dtype=float)
	distortion_array=np.zeros(len(current_array),dtype=float)
	output_dictionary={}

	# Storing the variable names
	size_array=['ss','sl','ls','ll']
	circuit_parameters={
		'len':0,
		'wid':0,
		'i_sin':1e-6,
		'v_1':1.0,
		'v_2':0.0,
		'v_b':0,
		'fund_1':1e9,
		'n_harm':15,
		'cir_temp':27
	}
	
	# Performing the analysis
	for resistor in resistor_dict:

		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)

		# Storing the resistor body voltage
		circuit_parameters['v_b']=resistor_dict[resistor]['v_body']

		for size in size_array:
			
			if size=='ss':
				circuit_parameters['wid']=resistor_dict[resistor]['w_min']
				circuit_parameters['len']=resistor_dict[resistor]['l_min']
			elif size=='sl':
				circuit_parameters['wid']=resistor_dict[resistor]['w_min']
				circuit_parameters['len']=resistor_dict[resistor]['l_max']
			elif size=='ls':
				circuit_parameters['wid']=resistor_dict[resistor]['w_max']
				circuit_parameters['len']=resistor_dict[resistor]['l_min']
			else:
				circuit_parameters['wid']=resistor_dict[resistor]['w_max']
				circuit_parameters['len']=resistor_dict[resistor]['l_max']

			# Creating the folder to store the outputs
			file_directory_current=file_directory_output+resistor+'/'+size+'/'
			if not os.path.exists(file_directory_current):
				os.makedirs(file_directory_current)
		
			# Opening the file
			filename_csv=file_directory_current+resistor+'_'+size+'_distortion.csv'
			f=open(filename_csv,'w')

			# Writing the first line in the csv file
			f.write('Current,Vout_fund_square,Vout_harmonics_square,Distortion,Vout_fund_dB,Vout_harmonics_dB,Distortion_dB\n')

			i=0
			for current in current_array:
			
				circuit_parameters['i_sin']=current
				resistance_dc,resistance_dict,distortion_dict,symmetry=write_extract(file_directory_netlist,circuit_parameters)
				
				f.write(str(current)+',')
				f.write(str(distortion_dict['vout_fund'])+',')
				f.write(str(distortion_dict['vout_harm'])+',')
				f.write(str(distortion_dict['distortion'])+',')
				f.write(str(distortion_dict['vout_fund_db'])+',')
				f.write(str(distortion_dict['vout_harm_db'])+',')
				f.write(str(distortion_dict['distortion_db'])+'\n')

				vout_fund_array[i]=distortion_dict['vout_fund_db']
				vout_harm_array[i]=distortion_dict['vout_harm_db']
				distortion_array[i]=distortion_dict['distortion_db']
				i+=1

			f.close()

			output_dictionary[size]={}
			output_dictionary[size]['fund']=vout_fund_array.copy()
			output_dictionary[size]['harm']=vout_harm_array.copy()
			output_dictionary[size]['distortion']=distortion_array.copy()
			
			# ---------- Plots ----------

			# Plot 1 - Vout_fund and Vout_extra
			figure()
			semilogx(current_array,vout_fund_array,color='green',label='Fundamental')
			semilogx(current_array,vout_harm_array,color='red',label='Harmonics')
			legend()
			grid()
			xlabel('Input Current')
			ylabel('Output Voltage ( in dB )')
			savefig(file_directory_current+resistor+'_'+size+'_vout.pdf')
			close()

			# Plot 2 - Distortion
			figure()
			semilogx(current_array,distortion_array,color='green',label='Distortion (dB)')
			legend()
			grid()
			xlabel('Input Current')
			ylabel('Distortion ( in dB )')
			savefig(file_directory_current+resistor+'_'+size+'_distortion.pdf')
			close()

		# ---------- Plots ----------
		file_directory_current=file_directory_output+resistor+'/'
		colour_dict={'ss':'green','sl':'red','ls':'blue','ll':'cyan'}

		# Plot 1 - Vout_fund and Vout_extra
		figure()
		for size in size_array:
			semilogx(current_array,output_dictionary[size]['fund'],color=colour_dict[size],linestyle='-',label='Fundamental '+size)
			semilogx(current_array,output_dictionary[size]['harm'],color=colour_dict[size],linestyle='--',label='Harmonics '+size)
		legend()
		grid()
		xlabel('Input Current')
		ylabel('Output Voltage ( in dB )')
		savefig(file_directory_current+resistor+'_vout.pdf')
		close()

		# Plot 2 - Distortion
		figure()
		for size in size_array:
			semilogx(current_array,output_dictionary[size]['distortion'],color=colour_dict[size],label=size)		
		legend()
		grid()
		xlabel('Input Current')
		ylabel('Distortion ( in dB )')
		savefig(file_directory_current+resistor+'_distortion.pdf')
		close()
		

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the distortion for all the resistors over extreme sizes and varying bias voltages but at a single current
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
# UPTO DATE
def MOS_Resistor_Distortion_single(file_directory_netlist,resistor_name,resistor_dict,file_directory_output):

	# Creating current array
	output_dictionary={}

	input_voltage_array=np.linspace(0,1,6)

	# Storing the variable names
	size_array=['ss','sl','ls','ll']
	circuit_parameters={
		'len':0,
		'wid':0,
		'i_sin':1e-4,
		'v_1':1.0,
		'v_2':0.0,
		'v_b':0,
		'fund_1':1e9,
		'n_harm':15,
		'cir_temp':27
	}
	
	# Writing the resistor name in the file
	write_resistor_name(file_directory_netlist,resistor_name)

	# Storing the resistor body voltage
	circuit_parameters['v_b']=resistor_dict['v_body']

	for size in size_array:
		print('\n\n\n\n\n\n\n Size:',size)
			
		if size=='ss':
			circuit_parameters['wid']=resistor_dict['w_min']
			circuit_parameters['len']=resistor_dict['l_min']
		elif size=='sl':
			circuit_parameters['wid']=resistor_dict['w_min']
			circuit_parameters['len']=resistor_dict['l_max']
		elif size=='ls':
			circuit_parameters['wid']=resistor_dict['w_max']
			circuit_parameters['len']=resistor_dict['l_min']
		else:
			circuit_parameters['wid']=resistor_dict['w_max']
			circuit_parameters['len']=resistor_dict['l_max']

		# Creating the folder to store the outputs
		file_directory_current=file_directory_output+resistor_name+'_'+size+'/'
		if not os.path.exists(file_directory_current):
			os.makedirs(file_directory_current)
		
		# Opening the file
		filename_csv=file_directory_current+resistor_name+'_'+size+'_distortion.csv'
		f=open(filename_csv,'w')

		# Writing the first line in the csv file
		f.write('v1,v2,Vout_fund_square,Vout_harmonics_square,Distortion,Vout_fund_dB,Vout_harmonics_dB,Distortion_dB\n')

		output_dictionary={}
		for v1 in input_voltage_array:
			output_dictionary[v1]={}
			v2_array=[]
			vout_fund_array=[]
			vout_harm_array=[]
			distortion_array=[]
			for v2 in input_voltage_array:
				if v1==v2:
					continue
				print('\n\n\n\n\n\n\n V1:',v1,' V2:',v2,' Start')
				circuit_parameters['v_1']=v1
				circuit_parameters['v_2']=v2
				resistance_dc,resistance_dict,distortion_dict,symmetry=write_extract(file_directory_netlist,circuit_parameters)
				
				f.write(str(v1)+',')
				f.write(str(v2)+',')
				f.write(str(distortion_dict['vout_fund'])+',')
				f.write(str(distortion_dict['vout_harm'])+',')
				f.write(str(distortion_dict['distortion'])+',')
				f.write(str(distortion_dict['vout_fund_db'])+',')
				f.write(str(distortion_dict['vout_harm_db'])+',')
				f.write(str(distortion_dict['distortion_db'])+'\n')

				v2_array.append(v2)
				vout_fund_array.append(distortion_dict['vout_fund_db'])
				vout_harm_array.append(distortion_dict['vout_harm_db'])
				distortion_array.append(distortion_dict['distortion_db'])
				print('\n\n\n\n\n\n\n V1:',v1,' V2:',v2,' End')
			
			output_dictionary[v1]['v2']=np.array(v2_array)
			output_dictionary[v1]['fund']=np.array(vout_fund_array)
			output_dictionary[v1]['harm']=np.array(vout_harm_array)
			output_dictionary[v1]['distortion']=np.array(distortion_array)

		f.close()

		# Plot 1
		colour_dict={0:'green',0.2:'red',0.4:'blue',0.6:'cyan',0.8:'lime',1.0:'maroon'}
		figure()
		for v1 in output_dictionary:
			plot(output_dictionary[v1]['v2'],output_dictionary[v1]['fund'],color=colour_dict[round(v1,2)],linestyle='-',label='Fundamental v1='+str(v1))
			plot(output_dictionary[v1]['v2'],output_dictionary[v1]['harm'],color=colour_dict[round(v1,2)],linestyle='-',label='Harmonics v1='+str(v1))
		legend()
		grid()
		xlabel('V2')
		ylabel('Output Voltage ( in dB )')
		savefig(file_directory_current+resistor_name+'_'+size+'_vout.pdf')
		close()

		# Plot 2
		figure()
		for v1 in output_dictionary:
			plot(output_dictionary[v1]['v2'],output_dictionary[v1]['distortion'],color=colour_dict[round(v1,2)],label='v1='+str(v1))
		legend()
		grid()
		xlabel('V2')
		ylabel('Distortion ( in dB )')
		savefig(file_directory_current+resistor_name+'_'+size+'_distortion.pdf')
		close()	
	


#---------------------------------------------------------------------------------------------------------------------------	
# Checking the symmetry
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Symmetry(file_directory_netlist,resistor_dict):

	# Storing the variables
	resistor_name=[]
	resistor_symmetry=[]
	
	# Performing the analysis
	for resistor in resistor_dict:
		
		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)
		
		# Choosing the width and length of the MOS Resistor
		wid=resistor_dict[resistor]['w_min']
		length=resistor_dict[resistor]['l_min']
		
		resistance_dc,resistance_ac,distortion,symmetry=write_extract(file_directory_netlist,length,wid,27,1e9)	# Finding resistance at 27o C
		resistor_name.append(resistor)
		resistor_symmetry.append(symmetry)
	
	print('\n\n\n\n\n Symmetry Analysis \n\n')
	for i in range(len(resistor_name)):
		print('\nResistor : ',resistor_name[i])
		print('Symmetry : ',resistor_symmetry[i])


#---------------------------------------------------------------------------------------------------------------------------	
# Performing DC Analysis
# Inputs: filenames for netlist files, resistor dictionary, output storing file directory
# Output: NONE
# UPTO DATE
def MOS_Resistor_DC_Analysis(file_directory_netlist,resistor_dict,file_directory_output):

	# Storing the variable names
	dc_input_array=np.linspace(0,1,5)
	size_array=['ss','sl','ls','ll']
	circuit_parameters={
		'len':0,
		'wid':0,
		'i_sin':1e-6,
		'v_1':0,
		'v_2':0,
		'v_b':0,
		'fund_1':1e9,
		'n_harm':15,
		'cir_temp':27
	}
	output_voltage_dictionary={}
	
	# Performing the analysis
	for resistor in resistor_dict:

		# Writing the resistor name in the file
		write_resistor_name(file_directory_netlist,resistor)
		print('\n\n Resistor : ', resistor)

		# Storing the resistor body voltage
		circuit_parameters['v_b']=resistor_dict[resistor]['v_body']
		
		# Choosing the width and length of the MOS Resistor
		for size in size_array:
			if size=='ss':
				circuit_parameters['wid']=resistor_dict[resistor]['w_min']
				circuit_parameters['len']=resistor_dict[resistor]['l_min']
			elif size=='sl':
				circuit_parameters['wid']=resistor_dict[resistor]['w_min']
				circuit_parameters['len']=resistor_dict[resistor]['l_max']
			elif size=='ls':
				circuit_parameters['wid']=resistor_dict[resistor]['w_max']
				circuit_parameters['len']=resistor_dict[resistor]['l_min']
			else:
				circuit_parameters['wid']=resistor_dict[resistor]['w_max']
				circuit_parameters['len']=resistor_dict[resistor]['l_max']

			# Creating the output path
			file_directory_current=file_directory_output+resistor+'/'+size+'/'
			print(file_directory_current)
			if not os.path.exists(file_directory_current):
				os.makedirs(file_directory_current)

			# Opening the file
			f=open(file_directory_current+resistor+'_'+size+'_'+'dc_resistance.csv','w')
			f.write('V1,V2,DC_Resistance\n')

			# Performing the analysis
			for v1 in dc_input_array:
				v2_array=[]
				resistance_array=[]
				for v2 in dc_input_array:
					if v1==v2:
						continue
					circuit_parameters['v_1']=v1
					circuit_parameters['v_2']=v2
					resistance_dc,resistance_ac_dict,distortion_dict,symmetry=write_extract(file_directory_netlist,circuit_parameters)
					f.write(str(v1)+',')
					f.write(str(v2)+',')
					f.write(str(resistance_dc)+'\n')
					
					v2_array.append(v2)
					resistance_array.append(resistance_dc)

				output_voltage_dictionary[v1]={}
				output_voltage_dictionary[v1]['v2_array']=np.array(v2_array)
				output_voltage_dictionary[v1]['resistance_array']=np.array(resistance_array)

			f.close()
			
			# Plotting the data
			figure()
			for v1 in output_voltage_dictionary:
				plot(output_voltage_dictionary[v1]['v2_array'],output_voltage_dictionary[v1]['resistance_array'],label='V1='+str(v1))
			xlabel('V2')
			ylabel('DC Resistance')
			grid()
			legend()
			savefig(file_directory_current+resistor+'_'+size+'_'+'resistance_plot.pdf')
			close()
			
"""



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


"""
# Code to do distortion analysis for multiple
write_directory_distortion='/home/ee18b028/Optimization/Simulation_Results/Resistance/Distortion_1_10v1/'
MOS_Resistor_Distortion(file_directory,resistor_dict_2,write_directory_distortion)
"""


"""
# Code to do distortion analysis for multiple
write_directory_distortion='/home/ee18b028/Optimization/Simulation_Results/Resistance/Distortion_1_10v1/'
MOS_Resistor_Distortion(file_directory,resistor_dict_2,write_directory_distortion)
"""

"""
# Code to do distortion analysis for single
write_directory_distortion2='/home/ee18b028/Optimization/Simulation_Results/Resistance/Distortion_Single_10v1/'
resistor_dict={'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':0.0}
MOS_Resistor_Distortion_single(file_directory,'rnodl_m',resistor_dict,write_directory_distortion2)
"""

"""
# Code to do DC Analysis
file_directory_output='/home/ee18b028/Optimization/Simulation_Results/Resistance/DC_Analysis_1_10v2/'
MOS_Resistor_DC_Analysis(file_directory,resistor_dict_2,file_directory_output)
"""

"""
# Code to frequency analysis
write_directory_fsweep='/home/ee18b028/Optimization/Simulation_Results/Resistance/FrequencySweep_2_10/'
MOS_Resistor_Frequency_Sweep(file_directory,resistor_dict_3,write_directory_fsweep)
"""

"""
# Code to do symmetry analysis
MOS_Resistor_Symmetry(file_directory,resistor_dict_2)
"""
