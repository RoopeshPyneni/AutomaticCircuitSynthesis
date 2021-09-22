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
# Extracting the DC from the file
# Inputs: Optimization_input_parameters
# Output: Dictionary with all the parameters

def extract_hb_param(filename):

	lines=extract_file(filename)
	
	frequency_array=[]
	vout_real_array=[]
	vout_img_array=[]
	flag=0

	for line in lines:
		if flag==1:
			char_real=(line.split()[1])[1:]
			char_img=(line.split()[2])[:-1]
			vout_real_array.append(float(char_real))
			vout_img_array.append(float(char_img))
			flag=0
		if 'freq' not in line:
			continue
		if 'sweep' in line.split()[1]:
			continue
		frequency_array.append(float(line.split()[1]))
		flag=1

	# Converting to numpy arrays
	frequency_array=np.array(frequency_array)
	vout_real_array=np.array(vout_real_array)
	vout_img_array=np.array(vout_img_array)

	vout_square_extra=0
	for i in range(len(frequency_array)):
		print('\n\n')
		print('i:',i)
		print('Frequency:',frequency_array[i])
		print('Vout Rea:',vout_real_array[i])
		print('Vout Img:',vout_img_array[i])
		if check_freq(frequency_array[i],0,1e6)==1:
			continue
		if check_freq(frequency_array[i],1e9,1e6)==1:
			resistance=1e6*np.sqrt(vout_real_array[i]**2+vout_img_array[i]**2)
			vout_square_normal=vout_real_array[i]**2+vout_img_array[i]**2
		else:
			vout_square_extra+=(vout_real_array[i]**2+vout_img_array[i]**2)
	distortion=vout_square_extra/(vout_square_extra+vout_square_normal)

	return resistance,distortion

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

def run_file():
	os.system('tcsh /home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh')	# This is the command to run the spectre file
	
#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def write_extract(filename_w,filename_e,filename_h,len,wid,temp):
	
	# Writing the tcsh file for Basic Analysis
	write_tcsh_file()

	# Writing to netlist file
	write_circuit_parameters(filename_w,len,wid,temp)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	resistance_dc=extract_dc_param(filename_e)

	# Extracting the HB Parameters
	resistance_ac,distortion=extract_hb_param(filename_h)
	
	return resistance_dc,resistance_ac,distortion

#===========================================================================================================================

"""
====================================================================================================================================================================================
------------------------------------------------------------ PLOTS -----------------------------------------------------------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will plot all the graphs
# Inputs  : Circuit_Parameters, Optimization_Input_Parameters
# Outputs : Extracted_Parameters

def plot_resistance(file_directory_plot,resistance_array,temp_array,len_array,wid_array,y_label):
	
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
		ylabel(y_label)
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
			loglog(len_array,resistance_array[i,:,k],label='Width = '+str(wid_array[k]))
		xlabel('Length')
		ylabel(y_label)
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
			loglog(wid_array,resistance_array[i,j,:],label='Length = '+str(len_array[j]))
		xlabel('Width')
		ylabel(y_label)
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
			loglog(len_array/wid_array[k],resistance_array[i,:,k],label='Width = '+str(wid_array[k]))
		xlabel('Length/Width')
		ylabel(y_label)
		grid()
		legend()
		savefig(file_directory+'/temp_'+str(temp_array[i])+'.pdf')
		close()



"""
====================================================================================================================================================================================
------------------------------------------------------------ TEMP CO FUNCTIONS -----------------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the slope and y-intercept
# Inputs: x and y coordinates of the points
# Output: slope, y-intercept

def calculate_slope(x,y):
	A = np.vstack([x, np.ones(len(x))]).T
	m, c = np.linalg.lstsq(A, y, rcond=None)[0]
	return m,c

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the temperature coefficient
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def temp_co_analysis(filename_w,filename_e,filename_h,resistor_list,file_directory):
	
	# Analysis performed to find the temperature coefficient
	temp_array=np.linspace(-40,120,17)
	resistance_array=np.zeros(17,dtype=float)
	filename_csv=file_directory+'/temp_coeff.csv'
	f=open(filename_csv,'w')

	# Writing the first line in the csv file
	f.write('Resistor Name,')
	for temp in temp_array:
		f.write(str(temp)+',')
	f.write('Temperature Coefficient\n')
	for resistor in resistor_list:
		f.write(resistor+',')
		if resistor=='rnwsti' or resistor=='rnwod':
			wid=5e-6
			len=10e-6
		else:
			wid=1e-6
			len=1e-6
		i=0
		for temp in temp_array:
			resistance_dc,resistance_ac,distortion=write_extract(filename_w,filename_e,filename_h,len,wid,temp)
			resistance_array[i]=resistance_dc
			i+=1
			f.write(str(resistance_dc)+',')
		temp_co,c=calculate_slope(temp_array,resistance_array)
		f.write(str(temp_co)+'\n')

	f.close()

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE DISTORTION -------------------------------------------------------------------------------------------
"""

def MOS_Resistor_Distortion():
	resistor='rpodwo'
	write_resistor_name(filename_w,resistor)
	resistance_dc,resistance_ac,distortion=write_extract(filename_w,filename_e,filename_h,1e-7,1e-7,27)

	print('Resistance DC:',resistance_dc)
	print('Resistance AC:',resistance_ac)
	print('Distortion   :',distortion)

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO SWEEP R vs W,L,T ----------------------------------------------------------------------------------------------
"""

def sweep_MOS_R():
	# Running the code 
	temp_array=np.linspace(-40,120,17)
	lw_start=60e-9
	lw_end=60e-7
	len_array=lw_start*np.logspace(0,np.log10(lw_end/lw_start),7)
	wid_array=len_array
	for resistor in resistor_list:
		print('Resistor name is : ',resistor)
		write_resistor_name(filename_w,resistor)

		file_directory_plot='/home/ee18b028/Optimization/Simulation_Results/Resistance/'+resistor
		if not os.path.exists(file_directory_plot):
			os.makedirs(file_directory_plot)

		l_temp=len(temp_array)
		l_len=len(len_array)
		l_wid=len(wid_array)

		resistance_dc_array=np.zeros((l_temp,l_len,l_wid),dtype=float)
		resistance_ac_array=np.zeros((l_temp,l_len,l_wid),dtype=float)
		distortion_array=np.zeros((l_temp,l_len,l_wid),dtype=float)

		for i in range(l_temp):
			for j in range(l_len):
				for k in range(l_wid):
					print('\n\ni=',i,'j=',j,'k=',k)
					resistance_dc_array[i,j,k],resistance_ac_array[i,j,k],distortion_array[i,j,k]=write_extract(filename_w,filename_e,filename_h,len_array[j],wid_array[k],temp_array[i])

		plot_resistance(file_directory_plot+'/DC_Resistance',resistance_dc_array,temp_array,len_array,wid_array,'DC Resistance')
		plot_resistance(file_directory_plot+'/AC_Resistance',resistance_ac_array,temp_array,len_array,wid_array,'AC Resistance')
		plot_resistance(file_directory_plot+'/Distortion',distortion_array,temp_array,len_array,wid_array,'Distortion')




"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""


# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/test/resistor_test_2'
filename_w=file_directory+'/circ.scs'
filename_e=file_directory+'/dc.out'
filename_h=file_directory+'/circ.raw/hb_test.fd.pss_hb'

# Creating the temperature, length, and width arrays
resistor_list=['rppolywo','rppolyl','rpodwo','rpodl','rnwsti','rnwod','rnpolywo','rnpolyl','rnodwo','rnodl']
#resistor_list=['rppolywo_m','rppolyl_m','rpodwo_m','rpodl_m','rnwsti_m','rnwod_m','rnpolywo_m','rnpolyl_m','rnodwo_m','rnodl_m']

write_directory='/home/ee18b028/Optimization/Simulation_Results/Resistance/Temp Coefficient'
temp_co_analysis(filename_w,filename_e,filename_h,resistor_list,write_directory)