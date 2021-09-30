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
			if 'R1a' in line:
				char_real=(line.split()[1])[1:]
				char_img=(line.split()[2])[:-1]
				vout_a_real_array.append(-1*float(char_img))
				vout_a_img_array.append(float(char_real))
			
			# Extracting R1b node values
			elif 'R1b' in line:
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

	# Calculating symmetry
	symmetry=0
	for i in range(len(frequency_array)):
		if check_freq(frequency_array[i],freq,freq/1000)==1:
			if vout_a_real_array[i]==-1*vout_b_real_array[i] and vout_a_img_array[i]==-1*vout_b_img_array[i]:
				symmetry=1
			break
	
	# Calculating AC Resistance and Impedance
	ac_resistance=0
	m_impedance=0
	p_impedance=0
	r_impedance=0
	i_impedance=0
	for i in range(len(frequency_array)):
		if check_freq(frequency_array[i],freq,freq/1000)==1:
			ac_resistance,m_impedance,p_impedance,r_impedance,i_impedance=calculate_impedance(i_cur,vout_a_real_array[i],vout_a_img_array[i],vout_b_real_array[i],vout_b_img_array[i])
			break
	
	resistance_dict={
		'AC_Resistance':ac_resistance,
		'Magnitude':m_impedance,
		'Phase':p_impedance,
		'Real':r_impedance,
		'Imaginary':i_impedance
	}

	return resistance_dict,distortion,symmetry

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

	# Calculating the AC Resistance
	Ca_Cb=-1*b_real/a_real
	a=a_real/cur
	b=a_img/cur
	c=a/(a**2+b**2)
	Resistance=(1+Ca_Cb)/c

	return Resistance,Z_mag,Z_ph,Z_real,Z_img



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
def write_circuit_parameters(filename,len,wid,temp,freq,i_cur):
	
	# Creating a dictionary of the parameters
	write_dict={'len':len,'wid':wid,'cir_temp':temp,'fund_1':freq,'i_sin':i_cur}

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
def write_extract(file_directory,length,wid,temp,freq,i_cur):

	# Getting the filenames
	filename_w=file_directory+'/circ.scs'
	filename_e=file_directory+'/dc.out'
	filename_h=file_directory+'/circ.raw/hb_test.fd.pss_hb'

	# Writing the tcsh file for Basic Analysis
	write_tcsh_file(file_directory)

	# Writing to netlist file
	write_circuit_parameters(filename_w,length,wid,temp,freq,i_cur)

	# Running netlist file
	run_file()

	# Extracting the Basic Parameters
	resistance_dc=extract_dc_param(filename_e)

	# Extracting the HB Parameters
	resistance_dict,distortion,symmetry=extract_hb_param(filename_h,freq,i_cur)
	
	return resistance_dc,resistance_dict,distortion,symmetry

#===========================================================================================================================


"""
====================================================================================================================================================================================
------------------------------------------------------------ FREQUENCY SWEEP ANALYSIS ----------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating AC Resistance as a function of frequency
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Frequency_Sweep1(file_directory_netlist,resistor_dict,file_directory):

	file_directory_ss=file_directory+'/ss'
	file_directory_sl=file_directory+'/sl'
	file_directory_ls=file_directory+'/ls'
	file_directory_ll=file_directory+'/ll'

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_ss):
		os.makedirs(file_directory_ss)

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_sl):
		os.makedirs(file_directory_sl)

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_ls):
		os.makedirs(file_directory_ls)

	# Creating the folder to store the outputs
	if not os.path.exists(file_directory_ll):
		os.makedirs(file_directory_ll)
	
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
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory_ss+'/FrequencyPlot_'+resistor+'.pdf')
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
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory_sl+'/FrequencyPlot_'+resistor+'.pdf')
		close()

		# ---------------------------- Case 3 --------------------------------
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
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory_ls+'/FrequencyPlot_'+resistor+'.pdf')
		close()

		# ---------------------------- Case 4 --------------------------------
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
		xlabel('Frequency')
		ylabel('Resistance')
		grid()
		legend()
		savefig(file_directory_ll+'/FrequencyPlot_'+resistor+'.pdf')
		close()

	f1.close()
	f2.close()
	f3.close()
	f4.close()

#---------------------------------------------------------------------------------------------------------------------------	
# Calculating AC Resistance as a function of frequency
# Inputs: filenames for netlist files, resistor list, output storing file directory
# Output: NONE
def MOS_Resistor_Frequency_Sweep(file_directory_netlist,resistor_dict,file_directory):

	size_types=['ss','sl','ls','ll']
	for size in size_types:
			
		# Getting the frequency array
		freq_array=np.logspace(8,10,11)

		# Arrays to store values
		resistance_ac_array=np.zeros(len(freq_array),dtype=float)
		impedance_m_array=np.zeros(len(freq_array),dtype=float)
		impedance_p_array=np.zeros(len(freq_array),dtype=float)
		impedance_r_array=np.zeros(len(freq_array),dtype=float)
		impedance_i_array=np.zeros(len(freq_array),dtype=float)
		resistance_dc_array=np.ones(len(freq_array),dtype=float)

		
		# Performing the analysis
		for resistor in resistor_dict:
			
			file_directory_current=file_directory+'/'+size+'/'+resistor

			# Creating the folder to store the outputs
			if not os.path.exists(file_directory_current):
				os.makedirs(file_directory_current)

			# Opening the file
			filename_csv=file_directory_current+'/frequency_sweep.csv'
			f=open(filename_csv,'w')

			# Writing the first line in the csv file
			f.write('Frequency,AC_Resistance,Impedance_Mag,Impedance_Phase,Impedance_Real,Impedance_Img,DC_Resistance\n')
			
			# Writing the resistor name in the file
			print('\n\n Resistor : ', resistor)
			write_resistor_name(file_directory_netlist,resistor)
		
			# Choosing the width and length
			if size=='ss':
				wid=resistor_dict[resistor]['w_min']
				length=resistor_dict[resistor]['l_min']
			elif size=='ls':
				wid=resistor_dict[resistor]['w_max']
				length=resistor_dict[resistor]['l_min']
			elif size=='sl':
				wid=resistor_dict[resistor]['w_min']
				length=resistor_dict[resistor]['l_max']
			else:
				wid=resistor_dict[resistor]['w_max']
				length=resistor_dict[resistor]['l_max']
			
			i=0
			for freq in freq_array:
				# Running spectre
				resistance_dc,resistance_dict,distortion,symmetry=write_extract(file_directory_netlist,length,wid,27,freq,1e-6)	# Finding resistance at 27o C

				# Storing the values in arrays
				resistance_dc_array[i]=resistance_dc
				resistance_ac_array[i]=resistance_dict['AC_Resistance']
				impedance_m_array[i]=resistance_dict['Magnitude']
				impedance_p_array[i]=resistance_dict['Phase']
				impedance_r_array[i]=resistance_dict['Real']
				impedance_i_array[i]=resistance_dict['Imaginary']

				# Writing the values
				f.write(str(freq)+',')
				f.write(str(resistance_dict['AC_Resistance'])+',')
				f.write(str(resistance_dict['Magnitude'])+',')
				f.write(str(resistance_dict['Phase'])+',')
				f.write(str(resistance_dict['Real'])+',')
				f.write(str(resistance_dict['Imaginary'])+',')
				f.write(str(resistance_dc)+'\n')
				i+=1
			
			f.close()

			# ---------- Plots ----------

			# Plot 1 : AC Resistance
			figure()
			semilogx(freq_array,resistance_ac_array,color='green',label='AC Resistance')
			semilogx(freq_array,resistance_dc_array,color='red',label='DC Resistance')
			xlabel('Frequency')
			ylabel('Resistance')
			grid()
			legend()
			savefig(file_directory_current+'/AC_Resistance.pdf')
			close()

			# Plot 2 : Impedance ( Mag and Phase )
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
			savefig(file_directory_current+'/Impedance_M_P.pdf')
			close()

			# Plot 3 : Impedance ( Real and Imaginary )
			figure()
			subplot(2,1,1)
			semilogx(freq_array,impedance_r_array,color='green',label='Impedance Real')
			ylabel('Impedance Real')
			grid()
			legend()
			subplot(2,1,2)
			semilogx(freq_array,impedance_i_array,color='red',label='Impedance Imaginary')
			xlabel('Frequency')
			ylabel('Impedance Imaginary')
			grid()
			legend()
			savefig(file_directory_current+'/Impedance_R_I.pdf')
			close()

		
	

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
file_directory='/home/ee18b028/cadence_project/test/resistor_test_4'

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

resistor_dict_3={
	'rnodl_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':2e-6,'w_max':10e-6},
	'rnods_m':{'l_min':0.4e-6,'l_max':100e-6,'w_min':0.4e-6,'w_max':2e-6}
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

#""
# Code to frequency analysis
write_directory_fsweep='/home/ee18b028/Optimization/Simulation_Results/Resistance/FrequencySweep_30_9'
MOS_Resistor_Frequency_Sweep(file_directory,resistor_dict_2,write_directory_fsweep)
#"""

"""
# Code to do symmetry analysis
MOS_Resistor_Symmetry(file_directory,resistor_dict_2)
"""

"""
====================================================================================================================================================================================
------------------------------------------------------------ TEMP CO FUNCTIONS -----------------------------------------------------------------------------------------------------

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
