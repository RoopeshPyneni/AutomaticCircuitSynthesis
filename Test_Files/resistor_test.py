#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will contain different analysis that is performed to characterize the resistors in TSMC 65nm model file
"""
#===========================================================================================================================
import numpy as np
import fileinput
import os
from matplotlib import pylab
from pylab import *

"""
====================================================================================================================================================================================
------------------------------------------------------------ EXTRACTION FUNCTION ---------------------------------------------------------------------------------------------------
"""

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

#---------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
# Inputs: file name
# Output: array of lines
def extract_file(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()
	return lines

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
	
	# Creating arrays to store the values
	frequency_array=[]
	vout_a_real_array=[]
	vout_a_img_array=[]
	vout_b_real_array=[]
	vout_b_img_array=[]
	
	# Extracting the data from the HB Raw File
	flag=0
	for line in lines:
		# If flag=1, it means that the "freq" line is already over and we need to extract the vout values.
		if flag==1:
			
			# Extracting R1a node values
			if 'Rn1' in line:
				char_real=(line.split()[1])[1:]
				char_img=(line.split()[2])[:-1]
				vout_a_real_array.append(-1*float(char_img))
				vout_a_img_array.append(float(char_real))
			
			# Extracting R1b node values
			elif 'Rn2' in line:
				char_real=(line.split()[1])[1:]
				char_img=(line.split()[2])[:-1]
				vout_b_real_array.append(-1*float(char_img))
				vout_b_img_array.append(float(char_real))
				flag=0
			else:
				continue
		if 'freq' not in line:
			continue
		if 'sweep' in line.split()[1]:
			continue
		frequency_array.append(float(line.split()[1]))
		flag=1

	# Converting to numpy arrays
	frequency_array=np.array(frequency_array)
	vout_a_real_array=np.array(vout_a_real_array)
	vout_a_img_array=np.array(vout_a_img_array)
	vout_b_real_array=np.array(vout_b_real_array)
	vout_b_img_array=np.array(vout_b_img_array)

	# Calculating distortion
	vout_square_extra=0
	for i in range(len(frequency_array)):
		if check_freq(frequency_array[i],0,freq/1000)==1:
			continue
		if check_freq(frequency_array[i],freq,freq/1000)==1:
			resistance=1e6*np.sqrt(vout_a_real_array[i]**2+vout_a_img_array[i]**2)
			vout_square_normal=vout_a_real_array[i]**2+vout_a_img_array[i]**2
		else:
			vout_square_extra+=(vout_a_real_array[i]**2+vout_a_img_array[i]**2)
	distortion=vout_square_extra/vout_square_normal
	distortion_dict={}
	distortion_dict['vout_fund']=vout_square_normal
	distortion_dict['vout_harm']=vout_square_extra
	distortion_dict['distortion']=distortion
	distortion_dict['vout_fund_db']=10*np.log10(vout_square_normal)
	distortion_dict['vout_harm_db']=10*np.log10(vout_square_extra)
	distortion_dict['distortion_db']=10*np.log10(distortion)

	# Calculating symmetry
	symmetry=0
	for i in range(len(frequency_array)):
		if check_freq(frequency_array[i],freq,freq/1000)==1:
			if vout_a_real_array[i]==-1*vout_b_real_array[i] and vout_a_img_array[i]==-1*vout_b_img_array[i]:
				symmetry=1
			break
	
	# Calculating Impedance
	m_impedance=0
	p_impedance=0
	r_impedance=0
	i_impedance=0
	for i in range(len(frequency_array)):
		if check_freq(frequency_array[i],freq,freq/1000)==1:
			m_impedance,p_impedance,r_impedance,i_impedance=calculate_impedance(i_cur,vout_a_real_array[i],vout_a_img_array[i],vout_b_real_array[i],vout_b_img_array[i])
			break
	
	resistance_dict={
		'Magnitude':m_impedance,
		'Phase':p_impedance,
		'Real':r_impedance,
		'Imaginary':i_impedance
	}

	return resistance_dict,distortion_dict,symmetry

#---------------------------------------------------------------------------------------------------------------------------	
# Extracting the Impedances from the voltages of node R1a and R1b
# Inputs: current, a_real, a_img, b_real, b_img 
# Output: resistance, (magintude, phase, real, and imaginary ) part of the impedance
def calculate_impedance(cur,a_real,a_img,b_real,b_img):
	
	# Calculating the Impedance
	Z_real=(a_real-b_real)/cur
	Z_img=(a_img-b_img)/cur
	Z_mag=np.sqrt(Z_real**2+Z_img**2)
	Z_ph=np.arctan(Z_img/Z_real)*57.29577

	return Z_mag,Z_ph,Z_real,Z_img



"""
====================================================================================================================================================================================
------------------------------------------------------------ FILE WRITE FUNCTIONS --------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that modifies the .scs file
# Inputs  : circuit_parameters, optimization input parameters
# Outputs : NONE
def write_resistor_name(file_directory,resistor_name):
	
	# We will write the new values to the Basic Circuit
	filename=file_directory+'/circ.scs'
	f=open(filename,'r+')
	s=''
	for line in fileinput.input(filename):
		if "R1" in line:	# Checking for a particular parameter in the .scs file
			line=line.split()
			if line[0]=='R1':
				line[4]=resistor_name
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
	filename='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
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
# Outputs : resistance_dc,resistance_dict,distortion_dict,symmetry
def write_extract(file_directory,circuit_parameters):

	# Getting the filenames
	filename_w=file_directory+'/circ.scs'
	filename_e=file_directory+'/dc.out'
	filename_h=file_directory+'/circ.raw/hb_test.fd.pss_hb'

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(file_directory)

	# Writing to netlist file
	write_circuit_parameters(filename_w,circuit_parameters)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	resistance_dc=extract_dc_param(filename_e)

	# Extracting the HB Parameters
	freq=circuit_parameters['fund_1']
	i_cur=circuit_parameters['i_sin']
	resistance_dict,distortion_dict,symmetry=extract_hb_param(filename_h,freq,i_cur)
	
	return resistance_dc,resistance_dict,distortion_dict,symmetry

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ FREQUENCY SWEEP ANALYSIS ----------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating AC Resistance as a function of frequency
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Frequency_Sweep(file_directory_netlist,resistor_dict,file_directory_output):

	# Storing the variable names
	size_types=['ss','sl','ls','ll']
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

	# Getting the frequency array
	freq_array=np.logspace(8,10,11)

	# Performing the analysis
	for resistor in resistor_dict:

		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)

		# Storing the resistor body voltage
		circuit_parameters['v_b']=resistor_dict[resistor]['v_body']

		for size in size_types:

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
			filename_csv=file_directory_current+resistor+'_'+size+'_frequency_sweep.csv'
			f=open(filename_csv,'w')

			# Writing the first line in the csv file
			f.write('Frequency,Impedance_Mag,Impedance_Phase,Impedance_Real,Impedance_Img,DC_Resistance\n')
			
			# Arrays to store values
			impedance_m_array=np.zeros(len(freq_array),dtype=float)
			impedance_p_array=np.zeros(len(freq_array),dtype=float)
			impedance_r_array=np.zeros(len(freq_array),dtype=float)
			impedance_i_array=np.zeros(len(freq_array),dtype=float)
			resistance_dc_array=np.ones(len(freq_array),dtype=float)
			
			i=0
			for freq in freq_array:
				
				circuit_parameters['fund_1']=freq

				# Running spectre
				resistance_dc,resistance_dict,distortion,symmetry=write_extract(file_directory_netlist,circuit_parameters)

				# Storing the values in arrays
				resistance_dc_array[i]=round(resistance_dc,2)
				impedance_m_array[i]=round(resistance_dict['Magnitude'],2)
				impedance_p_array[i]=round(resistance_dict['Phase'],2)
				impedance_r_array[i]=round(resistance_dict['Real'],2)
				impedance_i_array[i]=round(resistance_dict['Imaginary'],2)

				# Writing the values
				f.write(str(freq)+',')
				f.write(str(resistance_dict['Magnitude'])+',')
				f.write(str(resistance_dict['Phase'])+',')
				f.write(str(resistance_dict['Real'])+',')
				f.write(str(resistance_dict['Imaginary'])+',')
				f.write(str(resistance_dc)+'\n')
				i+=1
			
			f.close()

			# ---------- Plots ----------

			# Plot 1 : Impedance ( Mag and Phase )
			figure()
			subplot(2,1,1)
			semilogx(freq_array,impedance_m_array,color='green',label='Impedance Magnitude')
			semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
			legend()
			grid()
			ylabel('Impedance Magnitude')
			subplot(2,1,2)
			semilogx(freq_array,impedance_p_array,color='green',label='Impedance Phase')
			xlabel('Frequency')
			ylabel('Impedance Phase')
			grid()
			legend()
			savefig(file_directory_current+resistor+'_'+size+'_Impedance_M_P.pdf')
			close()

			# Plot 2 : Impedance ( Real and Imaginary )
			figure()
			subplot(2,1,1)
			semilogx(freq_array,impedance_r_array,color='green',label='Impedance Real')
			semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
			ylabel('Impedance Real')
			grid()
			legend()
			subplot(2,1,2)
			semilogx(freq_array,impedance_i_array,color='green',label='Impedance Imaginary')
			xlabel('Frequency')
			ylabel('Impedance Imaginary')
			grid()
			legend()
			savefig(file_directory_current+resistor+'_'+size+'_Impedance_R_I.pdf')
			close()

		
"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE DISTORTION -------------------------------------------------------------------------------------------
"""

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
	

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE SYMMETRY ---------------------------------------------------------------------------------------------
"""

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

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE DC RESISTANCE ----------------------------------------------------------------------------------------
"""

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
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/test/resistor_test_5'

# Creating the resistor dictionary
resistor_dict_1={
	'rnodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':0.0},
	'rnods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':2e-6,'v_body':0.0},
	'rnodwo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6,'v_body':00},
	'rnpolyl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':0.0},
	'rnpolys_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6,'v_body':0.0},
	'rnpolywo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6,'v_body':0.0},
	'rnwod_m':{'l_min':9e-6,'l_max':100e-6,'w_min':1.8e-6,'w_max':10e-6,'v_body':0.0},
	'rnwsti_m':{'l_min':9e-6,'l_max':100e-6,'w_min':1.8e-6,'w_max':10e-6,'v_body':0.0},
	'rpodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':1.0},
	'rpods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6,'v_body':1.0},
	'rpodwo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6,'v_body':1.0},
	'rppolyl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':1.0},
	'rppolys_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6,'v_body':1.0},
	'rppolywo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6,'v_body':1.0}
}

resistor_dict_2={
	'rnodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6,'v_body':0.0},
	'rnods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':2e-6,'v_body':0.0}
}

resistor_dict_3={
	'rppolywo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6,'v_body':1.0}
}


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

#"""
# Code to frequency analysis
write_directory_fsweep='/home/ee18b028/Optimization/Simulation_Results/Resistance/FrequencySweep_2_10/'
MOS_Resistor_Frequency_Sweep(file_directory,resistor_dict_3,write_directory_fsweep)
#"""

"""
# Code to do symmetry analysis
MOS_Resistor_Symmetry(file_directory,resistor_dict_2)
"""
