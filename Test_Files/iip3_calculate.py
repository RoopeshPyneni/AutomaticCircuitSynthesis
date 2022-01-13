#==========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will extract the value of iip3 from hb analysis file
"""

#==========================================================================================================================
import numpy as np


"""
===========================================================================================================================
------------------------------------- CHARACTER TO NUMBER -----------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------
# Changing the values extracted as a string to a floating point value 
def valueName_to_value(value_name):

	# Checking if the last character of array is a string
	if value_name[-1].isalpha()==0:
		val=float(value_name)
		return val
	
	# Checking if the last character of array is a string
	if (value_name[-1]=='G' and value_name[-2]=='E') or (value_name[-1]=='g' and value_name[-2]=='e'):
		val=float(value_name[:-3])*1e6
		return val
		
	# Extracting the numerical part of the number 
	val=float(value_name[:-1])
	
	# Extracting the character that denotes the units ( i.e, millt, micro, nano, etc)
	mult_name=value_name[-1]
	mult=1.0
	
	# Calculating the value of the unit
	if mult_name=='M' or mult_name=='m':
		mult=1e-3
	elif mult_name=='U' or mult_name=='u':
		mult=1e-6
	elif mult_name=='N' or mult_name=='n':
		mult=1e-9
	elif mult_name=='P' or mult_name=='p':
		mult=1e-12
	elif mult_name=='F' or mult_name=='f':
		mult=1e-15
	elif mult_name=='G' or mult_name=='g':
		mult=1e9
	else:
		mult=1.0
		
	val=val*mult
	return val
	
#--------------------------------------------------------------------------------------------------------------------------
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


"""
===========================================================================================================================
------------------------------------- IIP3 CALCULATION --------------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the IIP3 from a single point
def calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin):

	# Calculating values in log scale
	vout_fund_log=20*np.log10(vout_fund_mag)
	vout_im3_log=20*np.log10(vout_im3_mag)

	# Calculating iip3
	iip3=pin+(0.5*(vout_fund_log-vout_im3_log))

	return iip3,vout_fund_log,vout_im3_log,pin


"""
===========================================================================================================================
------------------------------------- FILE VALUE EXTRACTION ---------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting Vout magnitude of fundamental and im3 from file ( for hb_sweep )
def extract_vout_magnitude(file_name,fund_1,fund_2):

	lines=extract_file(file_name)
	
	# Getting f_error
	f_im3=2*fund_2-fund_1
	f_error=(fund_2-fund_1)/100


	flag=0
	flag_fun=0
	flag_im3=0
	flag_test=0

	while 1:
		if len(lines[0].split())<2:
			lines=lines[1:]
		
		elif 'freq' in lines[0].split()[0] and flag==0:
			flag=1
			lines=lines[1:]
		
		elif 'freq' in lines[0].split()[0] and flag==1:
			if flag_fun==0 and check_freq(float(lines[0].split()[1]),fund_2,f_error)==1 :
				
				#Extracting Vout for fundamental
				flag_test=1
				while flag_test==1:
					if 'Vout' in lines[0].split()[0]:
						flag_test=0
						vout_fund=extract_vout(lines[0])
					else:
						lines=lines[1:]
				flag_fun=1
			
			elif flag_im3==0 and check_freq(float(lines[0].split()[1]),f_im3,f_error)==1 :
				
				#Extracting Vout for im3
				flag_test=1
				while flag_test==1:
					if 'Vout' in lines[0].split()[0]:
						flag_test=0
						vout_im3=extract_vout(lines[0])
					else:
						lines=lines[1:]
				flag_im3=1
			lines=lines[1:]
			
			if flag_fun==1 and flag_im3==1:
				break
		
		else:
			lines=lines[1:]

	return vout_fund,vout_im3

#--------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
def extract_file(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()
	return lines

#--------------------------------------------------------------------------------------------------------------------------	
# Extracts Vout_magnitude from hb,pss file line
def extract_vout(lines):
	
	# Extracting Vout Magnitude
	lines=lines.split()
	char_r=lines[1].split('(')[1]
	char_i=lines[2].split(')')[0]

	# Converting string to floating point value
	vout_r=valueE_to_value(char_r)
	vout_i=valueE_to_value(char_i)
	
	# Calculating the magnitude of the output
	vout_mag=np.sqrt(vout_r*vout_r+vout_i*vout_i)

	return vout_mag

#--------------------------------------------------------------------------------------------------------------------------	
# Checks if the frequency is within range ( within (target-error,target+error) )
def check_freq(f_test,f_target,f_error):
	
	if f_test<f_target+f_error and f_test>f_target-f_error:
		return 1
	else:
		return 0


"""
===========================================================================================================================
------------------------------------- MAIN EXTRACTION FUNCTION ------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will extract the value of iip3 from the file
def extract_iip3(file_name,fund_1,fund_2,pin):
	
	# Extracting Vout Magnitude
	vout_fund_mag,vout_im3_mag=extract_vout_magnitude(file_name,fund_1,fund_2)

	# Calculating the iip3
	iip3_extracted_parameters={}
	iip3_extracted_parameters['iip3_dbm'],iip3_extracted_parameters['iip3_fund'],iip3_extracted_parameters['iip3_im3'],iip3_extracted_parameters['iip3_pin']=calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin)

	return iip3_extracted_parameters


"""
===========================================================================================================================
------------------------------------- MAIN PROGRAM ------------------------------------------------------------------------
"""

file_name='/home/ee18b028/Optimization/Simulation_Results/CS_LNA_Io/'+'iip3_hb'+'/circ.raw/hb_test.fd.qpss_hb'
fund_1=1000e6
fund_2=1001e6
pin=-65

iip3_parameters=extract_iip3(file_name,fund_1,fund_2,pin)

#==========================================================================================================================
