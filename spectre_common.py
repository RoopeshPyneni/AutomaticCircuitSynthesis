#==========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will contain the functions to write, run, and read from the spectre files for CS LNA
"""

#==========================================================================================================================
import numpy as np
import os
import common_functions as cf # type: ignore
import pandas as pd


"""
===========================================================================================================================
------------------------------------- OPTIMIZATION FUNCTIONS --------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This is the ramp function
def ramp_func(x):
	if x>0:
		return x
	else:
		return 0


"""
===========================================================================================================================
------------------------------------ MOSFET EXTRACTION --------------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that extracts the MOSFET File Parameeters
def calculate_mos_parameters(circuit_initialization_parameters):
	
	# Setting Lmin and Vdd
	Lmin=circuit_initialization_parameters['MOS']['Lmin']
	vdd=circuit_initialization_parameters['MOS']['Vdd']
	cox=circuit_initialization_parameters['MOS']['cox']
	un=circuit_initialization_parameters['MOS']['un']
	vt=circuit_initialization_parameters['MOS']['vt']

	# Extracting From File
	mos_parameters = {'un':un,'cox':cox,'vt':vt,'Lmin':Lmin,'vdd':vdd}
	
	# Printing the MOSFET Parameters
	cf.print_MOS_parameters(mos_parameters)

	return mos_parameters


"""
===========================================================================================================================
------------------------------------- EXTRACTION FUNCTION -----------------------------------------------------------------
"""

#==========================================================================================================================
#------------------------------------ Character to Real Number Functions --------------------------------------------------

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


#==========================================================================================================================
#------------------------------------ Other File Extraction Functions -----------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------
# Extracting the files as an array of lines
def extract_file(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()
	return lines


#==========================================================================================================================
#------------------------------------ Basic File Extraction Functions -----------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the gain and angle from the vout and vin values
def calculate_gain_phase(vout_re,vout_im,vin_re,vin_im):
	
	# Calculating gain_dB
	gain=(vout_re**2+vout_im**2)/(vin_re**2+vin_im**2)
	gain_db=10*np.log10(gain)

	# Calculating the phase of vout and vin
	if vout_re>0:
		vout_phase=np.arctan(vout_im/vout_re)*180/np.pi
	else:
		vout_phase=180+np.arctan(vout_im/vout_re)*180/np.pi
	if vin_re>0:
		vin_phase=np.arctan(vin_im/vin_re)*180/np.pi
	else:
		vin_phase=180+np.arctan(vin_im/vin_re)*180/np.pi
	
	# Calculating the phase of the gain
	phase=vout_phase-vin_phase
	while phase<-180:
		phase+=180
	while phase>180:
		phase-=180

	return gain_db,phase


#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the value of K from the SP 
def calculate_k(s11_db,s12_db,s21_db,s22_db,s11_ph,s12_ph,s21_ph,s22_ph):
	
	# Calculating the magnitude in normal scale
	s11_mag=10**(s11_db/20)
	s12_mag=10**(s12_db/20)
	s21_mag=10**(s21_db/20)
	s22_mag=10**(s22_db/20)

	# Calculating the values in a+ib format
	s11_real=s11_mag*np.cos(s11_ph)
	s12_real=s12_mag*np.cos(s12_ph)
	s21_real=s21_mag*np.cos(s21_ph)
	s22_real=s22_mag*np.cos(s22_ph)

	s11_img=s11_mag*np.sin(s11_ph)
	s12_img=s12_mag*np.sin(s12_ph)
	s21_img=s21_mag*np.sin(s21_ph)
	s22_img=s22_mag*np.sin(s22_ph)

	# Calculating delta squared
	delta_squared=(s11_real*s22_real-s12_real*s21_real+s12_img*s21_img-s11_img*s22_img)**2+(s11_real*s22_img+s22_real*s11_img-s12_real*s21_img-s21_real*s12_img)**2

	# Calculating k
	k=(1-(s11_mag**2)-(s22_mag**2)+delta_squared)/(2*s21_mag*s12_mag)

	return k

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the value of Zin_R and Zin_I from the SP 
def calculate_Z(s11_db,s11_ph):
	
	# Calculating the magnitude in normal scale
	s11_mag=10**(s11_db/20)

	# Calculating the values in a+ib format
	s11_real=s11_mag*np.cos(s11_ph)
	s11_img=s11_mag*np.sin(s11_ph)

	# Calculating the outputs
	denominator=(1-s11_real)**2+s11_img**2
	Zin_R=50*(1-s11_real**2-s11_img**2)/denominator
	Zin_I=50*(2*s11_img)/denominator

	return Zin_R,Zin_I


#==========================================================================================================================
#------------------------------------ IIP3 Extraction Functions -----------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the IIP3 from a single point
def calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin):

	# Calculating values in log scale
	vout_fund_log=20*np.log10(vout_fund_mag)
	vout_im3_log=20*np.log10(vout_im3_mag)

	# Calculating iip3
	iip3=pin+(0.5*(vout_fund_log-vout_im3_log))

	return iip3,vout_fund_log,vout_im3_log,pin

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the IIP3 after extraction of Vout data
def calculate_iip3_multiple_points(n_pin,n_points,vout_fund_mag,vout_im3_mag,pin):

	# Calculating values in log scale
	vout_fund_log=20*np.log10(vout_fund_mag)
	vout_im3_log=20*np.log10(vout_im3_mag)
	
	# Creating arrays for slopes and y-intercepts of fundamental and im3 components
	fund_slope=np.zeros(n_pin+1-n_points,dtype=float)
	fund_intercept=np.zeros(n_pin+1-n_points,dtype=float)
	im3_slope=np.zeros(n_pin+1-n_points,dtype=float)
	im3_intercept=np.zeros(n_pin+1-n_points,dtype=float)

	# Calculating the slopes and y-intercepts
	for i in range(n_pin+1-n_points):
		fund_slope[i],fund_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_fund_log[i:i+n_points-1])
		im3_slope[i],im3_intercept[i]=calculate_slope(pin[i:i+n_points-1],vout_im3_log[i:i+n_points-1])
	
	# Finding the best points for iip3 calculation
	best_point=calculate_best_iip3_point(fund_slope,im3_slope)
	
	# Calculating the iip3 given the slope and y-intercept of fundamental and im3
	iip3=(im3_intercept[best_point]-fund_intercept[best_point])/(fund_slope[best_point]-im3_slope[best_point])

	return iip3,im3_intercept[best_point],im3_slope[best_point],fund_intercept[best_point],fund_slope[best_point]

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the slope and y-intercept
def calculate_slope(x,y):
	A = np.vstack([x, np.ones(len(x))]).T
	m, c = np.linalg.lstsq(A, y, rcond=None)[0]
	return m,c

#--------------------------------------------------------------------------------------------------------------------------	
# Calculating the point with slope closest to 1dB/dB for fund and 3dB/dB for im3
def calculate_best_iip3_point(fund_slope,im3_slope):
	
	# Getting the length of the array
	l=len(fund_slope)

	# Choosing the best point as the first point
	best_point=0
	best_error=(fund_slope[0]-1)**2+(im3_slope[0]-3)**2

	for i in range(1,l):
		error=(fund_slope[i]-1)**2+(im3_slope[i]-3)**2
		if error<best_error:	# Checking if the current point is better than the best point
			best_point=i
			best_error=error
	return best_point

#--------------------------------------------------------------------------------------------------------------------------	
# Checks if the frequency is within range ( within (target-error,target+error) )
def check_freq(f_test,f_target,f_error):
	if f_test<f_target+f_error and f_test>f_target-f_error:
		return 1
	else:
		return 0

#==========================================================================================================================


"""
===========================================================================================================================
------------------------------------- FILE WRITE FUNCTIONS ----------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Command that returns the string that has to be printed in the .scs file
def print_param(param_var,val):
	return "parameters "+param_var+'='+str(val)+'\n'

#-----------------------------------------------------------------
# Function that modifies tcsh file
def write_tcsh_file(circuit_initialization_parameters,optimiztion_type):
	
	filename=circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']
	f=open(filename,'r+')
	s=''
	
	s='#tcsh\n'
	s=s+'source ~/.cshrc\n'
	
	if optimiztion_type=='basic':
		s=s+'cd '+circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'\n'
	else:
		s=s+'cd '+circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'\n'
	
	s=s+'spectre circ.scs =log circ_log.txt\n'
	#s=s+'spectre circ.scs \n'
	s=s+'exit'
	
	f.truncate(0)
	f.write(s)
	f.close()

#==========================================================================================================================


"""
===========================================================================================================================
------------------------------------- SPECTRE RUNNING FUNCTIONS -----------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------	
# This function will run the shell commands to run Spectre
def run_file(circuit_initialization_parameters):
	os.system('cd /home/ee18b028/cadence_project')
	os.system('tcsh '+circuit_initialization_parameters['simulation']['standard_parameters']['tcsh'])	# This is the command to run the spectre file

"""
===========================================================================================================================
------------------------------------- TSMC INDUCTOR FUNCTIONS -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Finding the location of L
def find_target_L(data,L,start,end):

	n1=int((start+end+1)/2)
	
	if data.at[n1,'L']>=L and data.at[n1-1,'L']<=L:
		return n1
	elif data.at[n1,'L']>=L and data.at[n1-1,'L']>=L:
		return find_target_L(data,L,start,n1-1)
	else:
		return find_target_L(data,L,n1,end)


#---------------------------------------------------------------------------------------------------------------------------
# Finding the best point
def find_TSMC_Inductor(Q,L):

	# Print Finding Q and L
	print('Finding the following inductor:')
	print('Q=',Q)
	print('L=',L)
	# Reading the data
	file_directory='/home/ee18b028/Optimization/Simulation_Results/Inductor/Sweep2/'
	#file_directory='C:/Users/roope/Studies/IIT/Prof Projects/Circuit_Synthesis/Extra_Codes/'
	data=pd.read_csv(file_directory+'inductor_sweep_1.csv')
	n_rows=data.shape[0]

	# Finding the location of 0.98L and 1.02L
	n1=find_target_L(data,L*0.98,0,n_rows-1)
	n2=find_target_L(data,L*1.02,0,n_rows-1)

	# Finding the loss of all the points from n1 and n2
	loss_iter=[]
	for i in range(n1,1+n2):
		Qi=data.at[i,'Q']
		Li=data.at[i,'L']
		loss=((Qi-Q)/Q)**2
		loss+=((Li-L)/L)**2
		loss_iter.append(loss)

	# Finding the location of the smallest value in loss iter
	min_value=min(loss_iter)
	min_index=loss_iter.index(min_value)
	min_index+=n1

	return min_index,data.at[min_index,'Width'],data.at[min_index,'Radius'],data.at[min_index,'N_Turns'],data.at[min_index,'gdis'],data.at[min_index,'spc']

#-----------------------------------------------------------------      
# Function that converts resistance to length and width
def get_TSMC_resistor(resistance):
	sheet_resistance=124.45
	W_min=0.4e-6
	dW=0.0691e-6
	width=W_min-dW
	length=width*resistance/sheet_resistance
	
	return length,W_min

#-----------------------------------------------------------------      
# Function that converts capacitance to length and width for MOS capacitor
def calculate_MOS_capacitor(cap):
	cox=17.25*1e-3
	w_check=np.sqrt(cap/cox)
	if w_check>2e-5:
		length=2e-5
		width=cap/(cox*length)
		return width,length
	if w_check<1.2e-7:
		width=1.2e-7
		length=cap/(cox*width)
		return width,length
	return w_check,w_check


#-----------------------------------------------------------------      
# Function that converts capacitance to length and width for MOS capacitor
def calculate_mimcap(cap_required,cap_1):
	
	return 1+int(cap_required/cap_1)


#==========================================================================================================================
