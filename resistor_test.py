#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Resistor Test:
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

def extract_hb_param(filename,freq):

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
		if check_freq(frequency_array[i],0,1e6)==1:
			continue
		if check_freq(frequency_array[i],freq,1e6)==1:
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
def write_circuit_parameters(filename,len,wid,temp,freq):
	
	# Creating a dictionary of the parameters
	write_dict={'len':len,'wid':wid,'cir_temp':temp,'fund_1':freq}

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
def write_tcsh_file(file_directory):
	filename='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	f=open(filename,'r+')
	s=''
	
	s='#tcsh\n'
	s=s+'source ~/.cshrc\n'
	s=s+'cd '+file_directory+'\n'
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

def write_extract(file_directory,length,wid,temp,freq):

	# Getting the filenames
	filename_w=file_directory+'/circ.scs'
	filename_e=file_directory+'/dc.out'
	filename_h=file_directory+'/circ.raw/hb_test.fd.pss_hb'

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(file_directory)

	# Writing to netlist file
	write_circuit_parameters(filename_w,length,wid,temp,freq)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	resistance_dc=extract_dc_param(filename_e)

	# Extracting the HB Parameters
	resistance_ac,distortion=extract_hb_param(filename_h,freq)
	
	return resistance_dc,resistance_ac,distortion

#===========================================================================================================================


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
def temp_co_analysis(file_directory_netlist,resistor_list,file_directory):
	
	# Creating the folder to store the outputs
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
			
	# Creating the arrays to store the results
	temp_array=np.linspace(-40,120,17)
	resistance_array=np.zeros(17,dtype=float)
	
	# Opening the file
	filename_csv=file_directory+'/temp_coeff.csv'
	f=open(filename_csv,'w')

	# Writing the first line in the csv file
	f.write('Resistor Name,')
	for temp in temp_array:
		f.write(str(temp)+',')
	f.write('Temperature Coefficient\n')
	
	# Performing the analysis
	for resistor in resistor_list:
		
		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)
		f.write(resistor+',')
		
		# Choosing the width and length of the MOS Resistor
		freq=1e9
		if resistor=='rnwsti' or resistor=='rnwod':
			wid=5e-6
			length=10e-6
		else:
			wid=1e-6
			length=1e-6
		
		# Starting the sweep for temperature
		i=0
		for temp in temp_array:
			print('Temperature : ',temp)
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,temp,freq)	# Extracting the values
			resistance_array[i]=resistance_dc	
			i+=1
			f.write(str(resistance_dc)+',') # Writing to csv file
		temp_co,c=calculate_slope(temp_array,resistance_array)	# Calculating the slope of resistance vs temperature
		resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
		temp_co/=resistance_dc	# Dividing the slope by resistance at 27o C
		f.write(str(temp_co)+'\n')

	f.close()

"""
====================================================================================================================================================================================
------------------------------------------------------------ FREQUENCY SWEEP ANALYSIS ----------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating AC Resistance as a function of frequency
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Frequency_Sweep1(file_directory_netlist,resistor_list,file_directory):

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
	
	# Opening the file
	filename_csv=file_directory+'/frequency_sweep.csv'
	f=open(filename_csv,'w')

	# Getting the frequency array
	freq_array=np.logspace(8,10,11)

	# Arrays to store values
	resistance_ac_array=np.zeros(len(freq_array),dtype=float)

	# Writing the first line in the csv file
	f.write('Resistor Name,')
	for freq in freq_array:
		f.write(str(freq)+',')
	f.write('DC Resistance\n')
	
	# Performing the analysis
	for resistor in resistor_list:
		
		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)
		f.write(resistor+',')
		
		# Choosing the width and length of the MOS Resistor
		if resistor=='rnwsti' or resistor=='rnwod' or resistor=='rnwsti_m' or resistor=='rnwod_m':
			wid=5e-6
			length=10e-6
		else:
			wid=3e-6
			length=3e-6
		
		i=0
		for freq in freq_array:
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
			resistance_ac_array[i]=resistance_ac
			i+=1
			f.write(str(resistance_ac)+',')
		f.write(str(resistance_dc)+'\n')

		resistance_dc_array=resistance_dc*np.ones(len(resistance_ac_array),dtype=float)
		figure()
		semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
		semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
		grid()
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/FrequencyPlot_'+resistor+'.pdf')
		close()

	f.close()

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating AC Resistance as a function of frequency
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Frequency_Sweep(file_directory_netlist,resistor_dict,file_directory):

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
	
	# Opening the file
	filename_csv1=file_directory+'/frequency_sweep_ss.csv'
	filename_csv2=file_directory+'/frequency_sweep_sl.csv'
	filename_csv3=file_directory+'/frequency_sweep_ls.csv'
	filename_csv4=file_directory+'/frequency_sweep_ll.csv'
	f1=open(filename_csv1,'w')
	f2=open(filename_csv2,'w')
	f3=open(filename_csv3,'w')
	f4=open(filename_csv4,'w')

	# Getting the frequency array
	freq_array=np.logspace(8,10,11)

	# Arrays to store values
	resistance_ac_array=np.zeros(len(freq_array),dtype=float)

	# Writing the first line in the csv file
	f1.write('Resistor Name,')
	f2.write('Resistor Name,')
	f3.write('Resistor Name,')
	f4.write('Resistor Name,')
	for freq in freq_array:
		f1.write(str(freq)+',')
		f2.write(str(freq)+',')
		f3.write(str(freq)+',')
		f4.write(str(freq)+',')
	f1.write('DC Resistance\n')
	f2.write('DC Resistance\n')
	f3.write('DC Resistance\n')
	f4.write('DC Resistance\n')
	
	# Performing the analysis
	for resistor in resistor_dict:
		
		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)
		f1.write(resistor+',')
		f2.write(resistor+',')
		f3.write(resistor+',')
		f4.write(resistor+',')
		
		# ---------------------------- Case 1 --------------------------------
		# Choosing the minimum width and minimum length 
		wid=resistor_dict[resistor]['w_min']
		length=resistor_dict[resistor]['l_min']
		i=0
		for freq in freq_array:
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
			resistance_ac_array[i]=resistance_ac
			i+=1
			f1.write(str(resistance_ac)+',')
		f1.write(str(resistance_dc)+'\n')

		resistance_dc_array=resistance_dc*np.ones(len(resistance_ac_array),dtype=float)
		figure()
		semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
		semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
		grid()
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/FrequencyPlot_ss_'+resistor+'.pdf')
		close()

		# ---------------------------- Case 2 --------------------------------
		# Choosing the minimum width and maximum length 
		wid=resistor_dict[resistor]['w_min']
		length=resistor_dict[resistor]['l_max']
		i=0
		for freq in freq_array:
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
			resistance_ac_array[i]=resistance_ac
			i+=1
			f2.write(str(resistance_ac)+',')
		f2.write(str(resistance_dc)+'\n')

		resistance_dc_array=resistance_dc*np.ones(len(resistance_ac_array),dtype=float)
		figure()
		semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
		semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
		grid()
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/FrequencyPlot_sl_'+resistor+'.pdf')
		close()

		# ---------------------------- Case 1 --------------------------------
		# Choosing the maximum width and minimum length 
		wid=resistor_dict[resistor]['w_max']
		length=resistor_dict[resistor]['l_min']
		i=0
		for freq in freq_array:
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
			resistance_ac_array[i]=resistance_ac
			i+=1
			f3.write(str(resistance_ac)+',')
		f3.write(str(resistance_dc)+'\n')

		resistance_dc_array=resistance_dc*np.ones(len(resistance_ac_array),dtype=float)
		figure()
		semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
		semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
		grid()
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/FrequencyPlot_ls_'+resistor+'.pdf')
		close()

		# ---------------------------- Case 1 --------------------------------
		# Choosing the minimum width and minimum length 
		wid=resistor_dict[resistor]['w_max']
		length=resistor_dict[resistor]['l_max']
		i=0
		for freq in freq_array:
			resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
			resistance_ac_array[i]=resistance_ac
			i+=1
			f4.write(str(resistance_ac)+',')
		f4.write(str(resistance_dc)+'\n')

		resistance_dc_array=resistance_dc*np.ones(len(resistance_ac_array),dtype=float)
		figure()
		semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
		semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
		grid()
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory+'/FrequencyPlot_ll_'+resistor+'.pdf')
		close()

	f1.close()
	f2.close()
	f3.close()
	f4.close()

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO FIND THE DISTORTION -------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating the distortion
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Distortion(file_directory_netlist,resistor_list,file_directory):

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
	
	# Opening the file
	filename_csv=file_directory+'/distortion.csv'
	f=open(filename_csv,'w')

	# Writing the first line in the csv file
	f.write('Resistor Name,')
	f.write('Distortion\n')
	
	# Performing the analysis
	for resistor in resistor_list:
		
		# Writing the resistor name in the file
		print('\n\n Resistor : ', resistor)
		write_resistor_name(file_directory_netlist,resistor)
		f.write(resistor+',')
		
		# Choosing the width and length of the MOS Resistor
		freq=1e9
		if resistor=='rnwsti' or resistor=='rnwod':
			wid=5e-6
			length=10e-6
		else:
			wid=1e-6
			length=1e-6
		
		resistance_dc,resistance_ac,distortion=write_extract(file_directory_netlist,length,wid,27,freq)	# Finding resistance at 27o C
		print('Resistance DC:',resistance_dc)
		print('Resistance AC:',resistance_ac)
		print('Distortion   :',distortion)
		f.write(str(distortion)+'\n')

	f.close()

"""
====================================================================================================================================================================================
------------------------------------------------------------ CODE TO SWEEP R vs W,L,T ----------------------------------------------------------------------------------------------
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

#---------------------------------------------------------------------------------------------------------------------------	
# Sweeping W,L,T and plotting the Resistance Values
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def sweep_MOS_R(file_directory_netlist,resistor_list,file_directory):
	
	# Running the code 
	freq=1e9
	temp_array=np.linspace(-40,120,17)
	lw_start=60e-9
	lw_end=60e-7
	len_array=lw_start*np.logspace(0,np.log10(lw_end/lw_start),7)
	wid_array=len_array
	for resistor in resistor_list:
		print('Resistor name is : ',resistor)
		write_resistor_name(file_directory_netlist,resistor)

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
					resistance_dc_array[i,j,k],resistance_ac_array[i,j,k],distortion_array[i,j,k]=write_extract(file_directory_netlist,
					len_array[j],wid_array[k],temp_array[i],freq)

		plot_resistance(file_directory_plot+'/DC_Resistance',resistance_dc_array,temp_array,len_array,wid_array,'DC Resistance')
		plot_resistance(file_directory_plot+'/AC_Resistance',resistance_ac_array,temp_array,len_array,wid_array,'AC Resistance')
		plot_resistance(file_directory_plot+'/Distortion',distortion_array,temp_array,len_array,wid_array,'Distortion')




"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""


# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/test/resistor_test_3'

# Creating the temperature, length, and width arrays
resistor_list1=['rppolywo','rppolyl','rpodwo','rpodl','rnwsti','rnwod','rnpolywo','rnpolyl','rnodwo','rnodl']
resistor_list2=['rppolywo_m','rppolyl_m','rpodwo_m','rpodl_m','rnwsti_m','rnwod_m','rnpolywo_m','rnpolyl_m','rnodwo_m','rnodl_m']
resistor_list3=['rnpolywo','rnodl','rnwsti']
resistor_list4=['rnpolywo_m','rnodl_m','rnwsti_m']

resistor_dict_2={
	'rnodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6},
	'rnods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':2e-6},
	'rnodwo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6},
	'rnpolyl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6},
	'rnpolys_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6},
	'rnpolywo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6},
	'rnwod_m':{'l_min':9e-6,'l_max':100e-6,'w_min':1.8e-6,'w_max':10e-6},
	'rnwsti_m':{'l_min':9e-6,'l_max':100e-6,'w_min':1.8e-6,'w_max':10e-6},
	'rpodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6},
	'rpods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6},
	'rpodwo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6},
	'rppolyl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6},
	'rppolys_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':1e-6,'w_max':2e-6},
	'rppolywo_m':{'l_min':0.8e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':10e-6}
}


"""
# Code to run temp co analysis
write_directory_temp='/home/ee18b028/Optimization/Simulation_Results/Resistance/Temp Coefficient'
temp_co_analysis(file_directory,resistor_list1,write_directory)
"""

"""
# Code to do distortion analysis
write_directory_distortion='/home/ee18b028/Optimization/Simulation_Results/Resistance/Distortion'
MOS_Resistor_Distortion(file_directory,resistor_list2,write_directory_distortion)
"""

#"""
# Code to frequency analysis
write_directory_fsweep='/home/ee18b028/Optimization/Simulation_Results/Resistance/FrequencySweep'
MOS_Resistor_Frequency_Sweep(file_directory,resistor_list2,write_directory_fsweep)
#"""
