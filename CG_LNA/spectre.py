#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will contain the functions to write, run, and read from the spectre files for CG LNA
"""
#===========================================================================================================================
import os
import shutil
import numpy as np
import fileinput
import multiprocessing as mp
import CG_LNA.pre_optimization as pr # type: ignore
import spectre_common as sp # type: ignore
import copy
import sys

"""
====================================================================================================================================================================================
------------------------------------ CIRCUIT CLASS ---------------------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Creating a class for the circuit 
class Circuit():

	#-----------------------------------------------------------------------------------------------
	#---------------------------- Initialization Functions -----------------------------------------

	# Initialization of the object
	def __init__(self,circuit_initialization_parameters):
		self.initial_circuit_parameters={}
		self.circuit_parameters={}
		self.extracted_parameters={}

		self.initial_circuit_initialization_parameters=copy.deepcopy(circuit_initialization_parameters)
		# Getting the circuit name
		if self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_resistor':
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_r'
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_r'
		
		elif self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_capacitor':
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_rc'
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_rc'

		elif self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_inductor':
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_rcl'
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_rcl'

		elif self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='ideal':
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters'
			self.initial_circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb'
		
		else:
			sys.exit()
		
		self.circuit_initialization_parameters=copy.deepcopy(self.initial_circuit_initialization_parameters)
		self.mos_parameters=sp.calculate_mos_parameters(self.circuit_initialization_parameters)

	
	#-----------------------------------------------------------------------------------------------
	#---------------------------- Circuit Parameter Functions --------------------------------------

	# Running the circuit
	def run_circuit(self):
		self.circuit_parameters=get_final_circuit_parameters(self.initial_circuit_parameters,self.circuit_initialization_parameters)
		self.extracted_parameters=write_extract(self.circuit_parameters,self.circuit_initialization_parameters)

	# Running multiple circuits 
	def run_circuit_multiple(self,initial_circuit_parameters_dict):
		circuit_parameters_dict={}
		for i in initial_circuit_parameters_dict:
			circuit_parameters_dict[i]=get_final_circuit_parameters(initial_circuit_parameters_dict[i].copy(),self.circuit_initialization_parameters).copy()
		extracted_parameters_dict=write_extract_multiple_circuits(circuit_parameters_dict,self.circuit_initialization_parameters)
		return extracted_parameters_dict

	# Updating the circuit parameters and running the circuit
	def update_circuit(self,initial_circuit_parameters):
		self.initial_circuit_parameters=initial_circuit_parameters
		self.circuit_parameters=get_final_circuit_parameters(self.initial_circuit_parameters,self.circuit_initialization_parameters)
		self.extracted_parameters=write_extract(self.circuit_parameters,self.circuit_initialization_parameters)
	
	# Updating the circuit parameters and not running the circuit
	def update_circuit_parameters_1(self,initial_circuit_parameters):
		self.initial_circuit_parameters=initial_circuit_parameters
		self.circuit_parameters=get_final_circuit_parameters(self.initial_circuit_parameters,self.circuit_initialization_parameters)

	# Getting the initial circuit parameters
	def get_initial_circuit_parameters(self):
		return self.initial_circuit_parameters.copy()

	# Getting the circuit parameters
	def get_circuit_parameters(self):
		return self.circuit_parameters.copy()
	
	# Getting the extracted parameters
	def get_extracted_parameters(self):
		return self.extracted_parameters.copy()
	
	# Update circuit
	def update_circuit_state(self,initial_circuit_parameters,circuit_parameters,extracted_circuit_parameters):
		self.initial_circuit_parameters=initial_circuit_parameters.copy()
		self.circuit_parameters=circuit_parameters.copy()
		self.extracted_circuit_parameters=extracted_circuit_parameters.copy()
	
	#-----------------------------------------------------------------------------------------------
	#---------------------------- Simulation Parameter Functions -----------------------------------

	# Updating the simulation parameters
	def update_simulation_parameters(self,simulation_parameters):
		
		self.circuit_initialization_parameters=copy.deepcopy(self.initial_circuit_initialization_parameters)

		if 'standard_parameters' in simulation_parameters:
			for param_name in simulation_parameters['standard_parameters']:
				self.circuit_initialization_parameters['simulation']['standard_parameters'][param_name]=simulation_parameters['standard_parameters'][param_name]
		
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['n_harm']=self.circuit_initialization_parameters['simulation']['standard_parameters']['n_harm']
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=self.circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']

	# Writing the simulation parameters
	def write_simulation_parameters(self):
		write_simulation_parameters(self.circuit_initialization_parameters)
	
	# Updating the circuit temperature
	def update_temp(self,temp):
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=temp
	
	# Resetting the circuit temperature
	def reset_temp(self):
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=self.circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']	


	#-----------------------------------------------------------------------------------------------
	#---------------------------- Optimization Functions -------------------------------------------

	#-----------------------------------------------------------------------------------------------
	# This function calculates the loss for Io Optimization
	# KWARGS : This is used to choose the temperature and process if necessary
	def calc_loss(self,output_conditions,loss_weights,**kwargs):
		
		# Adding the start of the parameter string
		tp_string=''
		if 'process' in kwargs:
			tp_string=str(kwargs['temp'])+'_'+str(kwargs['process'])+'_'
		elif 'temp' in kwargs:
			tp_string=str(kwargs['temp'])+'_'

		# Extracted Values
		gain=self.extracted_parameters[tp_string+'gain_db']
		iip3=self.extracted_parameters[tp_string+'iip3_dbm']
		s11=self.extracted_parameters[tp_string+'s11_db']
		nf=self.extracted_parameters[tp_string+'nf_db']
		Io=self.extracted_parameters[tp_string+'Io']
		
		# Reference Values
		gain_ref=output_conditions['gain_db']
		iip3_ref=output_conditions['iip3_dbm']
		s11_ref=output_conditions['s11_db']
		nf_ref=output_conditions['nf_db']
		
		# Defining the weights to calculate Loss
		A1=loss_weights['gain_db']	# Weight for gain
		A2=loss_weights['iip3_dbm']	# Weight for iip3
		A3=loss_weights['s11_db']	# Weight for s11
		A4=loss_weights['nf_db']	# Weight for nf
		A5=loss_weights['Io']	# Weight for Io
		
		# Calculating Loss
		loss_gain=A1*sp.ramp_func(gain_ref-gain)
		loss_iip3=A2*sp.ramp_func(iip3_ref-iip3)
		loss_s11=A3*sp.ramp_func(s11-s11_ref)
		loss_nf=A4*sp.ramp_func(nf-nf_ref)
		loss_Io=A5*Io

		loss=loss_gain+loss_iip3+loss_s11+loss_nf+loss_Io
		loss_dict={
			'loss':loss,
			'loss_gain':loss_gain,
			'loss_iip3':loss_iip3,
			'loss_s11':loss_s11,
			'loss_nf':loss_nf,
			'loss_Io':loss_Io
		}
		
		return loss_dict
	
	#-----------------------------------------------------------------------------------------------
	# This function updates the values of circuit parameters by trying to minimize loss
	def update_circuit_parameters(self,circuit_parameters_slope,optimization_input_parameters,run_number,loss_iter,loss_type):
		
		# Getting the names of loss parameters that should be considered for circuit parameter updation
		change_loss_parameters=[]

		# Loss Type = 0
		if loss_type==0:
			if loss_iter['loss']==loss_iter['loss_Io']:
				change_loss_parameters=['loss_Io']
			else:
				change_loss_parameters=['loss_s11','loss_gain','loss_iip3','loss_nf']
		
		# Loss Type = 1
		elif loss_type==1:
			change_loss_parameters=['loss']

		# Loss Type = 2
		elif loss_type==2:
			for param in loss_iter:
				if param=='loss':
					continue
				if loss_iter[param]==0:
					continue
				change_loss_parameters.append(param)

		
		alpha=optimization_input_parameters['optimization'][run_number]['alpha']['value']

		# Calculating the value to update each parameter with
		for param_name in circuit_parameters_slope:
			
			# Calculating the Increment Value
			change=0
			for loss_name in change_loss_parameters:
				change+=circuit_parameters_slope[param_name][loss_name]
			change=change*(self.initial_circuit_parameters[param_name]**2)*alpha
		
		
			# Checking if the parameter is updated by a large value
			change_limit=0.25 # If the incremented value is more than +- change_limit*parameter_name, then we will limit the change
			if change>change_limit*self.initial_circuit_parameters[param_name]:
				change=change_limit*self.initial_circuit_parameters[param_name]
			if change<-1*change_limit*self.initial_circuit_parameters[param_name]:
				change=-1*change_limit*self.initial_circuit_parameters[param_name]
			
			# Updating circuit_parameters
			self.initial_circuit_parameters[param_name]=self.initial_circuit_parameters[param_name]-change
		
		self.update_circuit_parameters_1(self.initial_circuit_parameters)
			
	#---------------------------------------------------------------------------------------------------------------------------
	# Function to check the best solution
	def check_best_solution(self,optimization_results,loss_max):

		# Defining some values
		n_iter=optimization_results['n_iter']
		iter_min=-1

		# These arrays define which variable must be zero and which should be optimized
		zero_loss_array=['loss_s11','loss_gain','loss_iip3','loss_nf']
		minimize_loss_array=[]

		# Checking the first iter point
		total_loss=optimization_results['loss_iter'][-1]['loss']
		loss_Io_min=sum([optimization_results['loss_iter'][-1][key] for key in minimize_loss_array])
		if sum([optimization_results['loss_iter'][-1][key] for key in zero_loss_array])>loss_max:
			flag=0
		else:
			flag=1 # This means that we got a correct point
		
		# Checking other iterations
		for i in range(0,n_iter):
			if sum([optimization_results['loss_iter'][i][key] for key in zero_loss_array])>loss_max:
				if flag==1:
					continue
				if optimization_results['loss_iter'][i]['loss']<total_loss:
					iter_min=i
					total_loss=optimization_results['loss_iter'][i]['loss']
					loss_Io_min=sum([optimization_results['loss_iter'][i][key] for key in minimize_loss_array])
				continue

			if flag==0 or (flag==1 and sum([optimization_results['loss_iter'][i][key] for key in minimize_loss_array])<loss_Io_min):
				iter_min=i
				loss_Io_min=sum([optimization_results['loss_iter'][i][key] for key in minimize_loss_array])
				flag=1

		# Creating output dictionary
		opt_dict={}
		opt_dict['loss_max']=loss_max
		if flag==1:
			opt_dict['perfect_point']='Yes'
		else:
			opt_dict['perfect_point']='No'
		opt_dict['iter_number']=iter_min+1
		opt_dict['Io_loss']=loss_Io_min
		
		return opt_dict
	
	#-----------------------------------------------------------------------------------------------
	# This function will check which process and temperature is the worst case
	def check_worst_case(self,output_conditions,loss_weights):

		# Getting the temperature and process list
		temp_list=self.circuit_initialization_parameters['simulation']['standard_parameters']['temp_list']
		process_list=self.circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']

		# Getting the middle temperature and process
		temp_len=len(temp_list)
		process_len=len(process_list)

		temp_middle=temp_list[(temp_len-1)//2]
		process_middle=process_list[(process_len-1)//2]

		# Getting the list of parameters that need to be optimized
		zero_loss_array=['loss_s11','loss_gain','loss_iip3','loss_nf']

		# Choosing the centre point as the best point
		current_loss=self.calc_loss(output_conditions,loss_weights,temp=temp_middle,process=process_middle).copy()
		current_loss_total=sum([current_loss[param] for param in zero_loss_array])
		current_temp=temp_middle
		current_process=process_middle
		
		for temp in temp_list:
			for process in process_list:
				if temp==temp_middle and process==process_middle:
					continue
				loss=self.calc_loss(output_conditions,loss_weights,temp=temp,process=process)
				loss_total=sum([loss[param] for param in zero_loss_array])
				
				if loss_total>current_loss_total:
					current_loss_total=loss_total
					current_loss=loss.copy()
					current_temp=temp
					current_process=process
		
		return current_loss,current_temp,current_process

	#---------------------------------------------------------------------------------------------------------------------------
	# Function to perform pre optimization
	def pre_optimization(self,optimization_input_parameters,timing_results):
		pr.pre_optimization(self,optimization_input_parameters,timing_results)

	#-----------------------------------------------------------------------------------------------
	#---------------------------- Other Functions --------------------------------------------------

	#-----------------------------------------------------------------------------------------------
	# Calculating the iip3
	def calculate_iip3(self,n_pin,n_points,vout_fund_mag,vout_im3_mag,pin):
		return sp.calculate_iip3_multiple_points(n_pin,n_points,vout_fund_mag,vout_im3_mag,pin)


"""
===========================================================================================================================
------------------------------------- CONSTRAINTS FUNCTION ----------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------
# Getting the final circuit parameters
def get_final_circuit_parameters(initial_circuit_parameters,circuit_initialization_parameters):
	
	circuit_parameters=initial_circuit_parameters.copy()

	# ~~~~~~~~~~~~~~~ CONSTRAINTS CHECK ~~~~~~~~~~~~~~~

	
	# ~~~~~~~~~~~~~~~ GETTING EXTRA PARAMETERS ~~~~~~~~~~~~~~~

	# Getting the circuit type
	circuit_type=circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']

	# Calculating the number of fingers for the MOSFETs
	circuit_parameters['n_finger']=int(circuit_parameters['W']/circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max'])+1

	# Getting the real resistor parameters
	if circuit_type=='mos_resistor' or circuit_type=='mos_capacitor' or circuit_type=='mos_inductor':
		circuit_parameters['Resb_L'],circuit_parameters['Resb_W']=sp.get_TSMC_resistor(circuit_parameters['Rb'])
		circuit_parameters['Resd_L'],circuit_parameters['Resd_W']=sp.get_TSMC_resistor(circuit_parameters['Rd'])
		circuit_parameters['Resbias_L'],circuit_parameters['Resbias_W']=sp.get_TSMC_resistor(circuit_parameters['Rbias'])

	# Getting the real capacitors parameters
	if circuit_type=='mos_capacitor' or circuit_type=='mos_inductor':
		circuit_parameters['mf_cap1']=sp.calculate_mimcap(circuit_parameters['C1'],1e-11)	
		circuit_parameters['mf_cap2']=sp.calculate_mimcap(circuit_parameters['C2'],1e-11)	

		circuit_parameters['wid_cap2'],circuit_parameters['len_cap2']=sp.calculate_MOS_capacitor(circuit_parameters['C2'])
		
	# ~~~~~~~~~~~~~~~ CHANGING THE NAME ~~~~~~~~~~~~~~~

	cir_writing_dict={	
		'wid':'W',
		'cur0':'Io',
		'Resb':'Rb',
		'Resd':'Rd',
		'cap1':'C1',
		'cap2':'C2',
		'Resbias':'Rbias'
	}

	for param_name in cir_writing_dict:
		circuit_parameters[param_name]=circuit_parameters[cir_writing_dict[param_name]]
		del circuit_parameters[cir_writing_dict[param_name]]
	
	return circuit_parameters


"""
===========================================================================================================================
------------------------------------ EXTRACTION FUNCTION ------------------------------------------------------------------
"""

#==========================================================================================================================
#----------------------------------- Basic File Extraction Functions ------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the DC from the file
def extract_dc_param(circuit_initialization_parameters):

	# Getting the filename
	file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/dc.out'
	lines=sp.extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values from the required line
	extracted_parameters['vd']=sp.valueE_to_value(lines[1])
	extracted_parameters['vg']=sp.valueE_to_value(lines[2])
	extracted_parameters['vs']=sp.valueE_to_value(lines[3])

	extracted_parameters['i_source']=np.absolute(sp.valueE_to_value(lines[4]))
	extracted_parameters['v_source']=np.absolute(sp.valueE_to_value(lines[5]))
	extracted_parameters['p_source']=extracted_parameters['i_source']*extracted_parameters['v_source']

	extracted_parameters['Io']=sp.valueE_to_value(lines[6])
	extracted_parameters['gm1']=sp.valueE_to_value(lines[7])
	extracted_parameters['gds1']=sp.valueE_to_value(lines[8])
	extracted_parameters['vt']=sp.valueE_to_value(lines[9])
	extracted_parameters['vdsat']=sp.valueE_to_value(lines[10])

	extracted_parameters['cgd1']=np.absolute(sp.valueE_to_value(lines[11]))
	extracted_parameters['cgs1']=np.absolute(sp.valueE_to_value(lines[12]))

	return extracted_parameters

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the AC from the file
def extract_ac_param(circuit_initialization_parameters):

	# Getting the filename
	file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/ac.out'
	lines=sp.extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the values frim the required line
	extracted_parameters['freq']=sp.valueE_to_value(lines[0])
	vout_re=sp.valueE_to_value(lines[1])
	vout_im=sp.valueE_to_value(lines[2])
	vin_re=sp.valueE_to_value(lines[3])
	vin_im=sp.valueE_to_value(lines[4])
	extracted_parameters['gain_db'],extracted_parameters['gain_phase']=sp.calculate_gain_phase(vout_re,vout_im,vin_re,vin_im)

	return extracted_parameters

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the SP from the file
def extract_sp_param(circuit_initialization_parameters):

	# Getting the filename
	file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/sp.out'
	lines=sp.extract_file(file_name)
	
	extracted_parameters={}
	
	# Skipping the first few lines
	while 1:
		if 'format freq' not in lines[0]:
			lines=lines[1:]
		else:
			break
	lines=lines[4:]
	
	line1=lines[0].split()
	line2=lines[1].split()

	# Extracting the value of s-parameters in dB from the required line
	num_char_s11=line1[1].split(',')[0]
	num_char_s21=line1[3].split(',')[0]
	num_char_s12=line2[0].split(',')[0]
	num_char_s22=line2[2].split(',')[0]

	# Extracting the phase of the s-parameters
	num_char_s11_rad=float(line1[2])
	num_char_s21_rad=float(line1[4])
	num_char_s12_rad=float(line2[1])
	num_char_s22_rad=float(line2[3])

	# Calculating the final values
	extracted_parameters['s11_db']=sp.valueE_to_value(num_char_s11)
	extracted_parameters['s12_db']=sp.valueE_to_value(num_char_s12)
	extracted_parameters['s21_db']=sp.valueE_to_value(num_char_s21)
	extracted_parameters['s22_db']=sp.valueE_to_value(num_char_s22)

	extracted_parameters['k']=sp.calculate_k(extracted_parameters['s11_db'],extracted_parameters['s12_db'],extracted_parameters['s21_db'],extracted_parameters['s22_db'],
	num_char_s11_rad,num_char_s12_rad,num_char_s21_rad,num_char_s22_rad)

	extracted_parameters['Zin_R'],extracted_parameters['Zin_I']=sp.calculate_Z(extracted_parameters['s11_db'],num_char_s11_rad)

	return extracted_parameters

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting the Noise from the file
def extract_noise_param(circuit_initialization_parameters):

	# Getting the filename
	file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/noise.out'
	lines=sp.extract_file(file_name)
	
	extracted_parameters={}

	# Skipping the first few lines
	lines=lines[7:]
	lines=lines[0].split()

	# Extracting the value from the required line
	extracted_parameters['nf_db']=sp.valueE_to_value(lines[1])

	return extracted_parameters

#--------------------------------------------------------------------------------------------------------------------------
# Extracting all the output parameters from chi file
def extract_basic_parameters(circuit_initialization_parameters):
	
	# Extracting the outputs 
	extracted_parameters_dc=extract_dc_param(circuit_initialization_parameters)
	extracted_parameters_ac=extract_ac_param(circuit_initialization_parameters)
	extracted_parameters_sp=extract_sp_param(circuit_initialization_parameters)
	extracted_parameters_noise=extract_noise_param(circuit_initialization_parameters)
	
	# Storing the outputs in a single dictionary
	extracted_parameters={}
	
	for param_name in extracted_parameters_dc:
		extracted_parameters[param_name]=extracted_parameters_dc[param_name]
	for param_name in extracted_parameters_ac:
		extracted_parameters[param_name]=extracted_parameters_ac[param_name]
	for param_name in extracted_parameters_sp:
		extracted_parameters[param_name]=extracted_parameters_sp[param_name]
	for param_name in extracted_parameters_noise:
		extracted_parameters[param_name]=extracted_parameters_noise[param_name]
	
	return extracted_parameters


#==========================================================================================================================
#----------------------------------- IIP3 Extraction Functions ------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------	
# Extracting Vout magnitude of fundamental and im3 from file ( for hb_sweep )
def extract_vout_magnitude(file_name,circuit_initialization_parameters):

	lines=sp.extract_file(file_name)

	fund_1=circuit_initialization_parameters['simulation']['netlist_parameters']['fund_1']
	fund_2=circuit_initialization_parameters['simulation']['netlist_parameters']['fund_2']

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
			if flag_fun==0 and sp.check_freq(float(lines[0].split()[1]),fund_2,f_error)==1 :
				
				#Extracting Vout for fundamental
				flag_test=1
				while flag_test==1:
					if 'Vout' in lines[0].split()[0]:
						flag_test=0
						vout_fund=extract_vout(lines[0])
					else:
						lines=lines[1:]
				flag_fun=1
			
			elif flag_im3==0 and sp.check_freq(float(lines[0].split()[1]),f_im3,f_error)==1 :
				
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
# Extracts Vout_magnitude from hb,pss file line
def extract_vout(lines):
	
	# Extracting Vout Magnitude
	lines=lines.split()
	char_r=lines[1].split('(')[1]
	char_i=lines[2].split(')')[0]

	# Converting string to floating point value
	vout_r=sp.valueE_to_value(char_r)
	vout_i=sp.valueE_to_value(char_i)
	
	# Calculating the magnitude of the output
	vout_mag=np.sqrt(vout_r*vout_r+vout_i*vout_i)

	return vout_mag

#==========================================================================================================================



"""
===========================================================================================================================
------------------------------------ FILE WRITE FUNCTIONS -----------------------------------------------------------------
"""

#-----------------------------------------------------------------
# Function that modifies the .scs file
def write_circuit_parameters(circuit_parameters,circuit_initialization_parameters):
	
	# Getting the filenames
	filename1=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
	filename2=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.scs'

	# We will write the new values to the Basic Circuit
	f=open(filename1,'r+')
	s=''
	for line in fileinput.input(filename1):
		for param_name in circuit_parameters:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,sp.print_param(param_name,circuit_parameters[param_name]))	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

	# We will write the new values to the IIP3 Circuit
	f=open(filename2,'r+')
	s=''
	for line in fileinput.input(filename2):
		for param_name in circuit_parameters:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,sp.print_param(param_name,circuit_parameters[param_name]))	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that adds MOSFET Parameters to the netlist
def write_MOS_parameters(circuit_initialization_parameters):
	
	# write_dict will contain the mosfet values
	write_dict={
		'len':circuit_initialization_parameters['MOS']['Lmin'],
		'v_dd':circuit_initialization_parameters['MOS']['Vdd'],
	}
	process_corner=circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']

	# Getting the filenames
	filename1=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
	filename2=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.scs'

	# Writing the MOS Parameters to Basic File
	f=open(filename1,'r+')
	s=''
	write_check=1
	include_check=0

	# Replacing the lines of .scs file
	for line in fileinput.input(filename1):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+circuit_initialization_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,sp.print_param(param_name,write_dict[param_name]))

		if write_check==1:
			s=s+line

	f.truncate(0)
	f.write(s)
	f.close()

	# Writing the MOS Parameters to IIP3 File
	f=open(filename2,'r+')
	s=''
	write_check=1
	include_check=0

	# Replacing the lines of .scs file
	for line in fileinput.input(filename2):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+circuit_initialization_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,sp.print_param(param_name,write_dict[param_name]))
		
		if 'hb_test' in line and 'errpreset=conservative' in line and circuit_initialization_parameters['simulation']['standard_parameters']['conservative']=='NO':
			line_split=line.split()
			line=''
			for word in line_split:
				if 'errpreset=conservative' not in word:
					line=line+word+' '
			line=line+'\n'
		
		elif 'hb_test' in line and 'errpreset=conservative' not in line and circuit_initialization_parameters['simulation']['standard_parameters']['conservative']=='YES':
			line=line[:-1]+' errpreset=conservative \n'

		if write_check==1:
			s=s+line
	

	f.truncate(0)
	f.write(s)
	f.close()

#-----------------------------------------------------------------
# Function that adds Simulation Parameters
def write_simulation_parameters(circuit_initialization_parameters,circuit_name):
	
	# Adding simulation_parameters to write_dict
	write_dict={}
	for param_name in circuit_initialization_parameters['simulation']['netlist_parameters']:
		write_dict[param_name]=circuit_initialization_parameters['simulation']['netlist_parameters'][param_name]
	process_corner=circuit_initialization_parameters['simulation']['netlist_parameters']['process_corner']

	write_dict['len']=circuit_initialization_parameters['MOS']['Lmin']
	write_dict['v_dd']=circuit_initialization_parameters['MOS']['Vdd']

	# Getting the filenames
	if circuit_name=='basic':
		filename=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
	else:
		filename=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.scs'

	# Writing the simulation parameters to Basic File
	f=open(filename,'r+')
	s=''
	write_check=1
	include_check=0
	
	# Replacing the lines of .scs file
	for line in fileinput.input(filename):
		if "include " in line:	# This line is used to include the MOS file in the .scs file
			include_check=1
			write_check=0

		elif "include" not in line and include_check==1:
			s=s+circuit_initialization_parameters['MOS']['filename'][process_corner]
			include_check=0
			write_check=1
		
		for param_name in write_dict:	# This line is used to replace the MOS parameters and simulation_parameters
			if "parameters "+param_name+'=' in line:
				line=line.replace(line,sp.print_param(param_name,write_dict[param_name]))
		
		if circuit_name=='iip3':
			if 'hb_test' in line and 'errpreset=conservative' in line and circuit_initialization_parameters['simulation']['standard_parameters']['conservative']=='NO':
				line_split=line.split()
				line=''
				for word in line_split:
					if 'errpreset=conservative' not in word:
						line=line+word+' '
				line=line+'\n'
			
			elif 'hb_test' in line and 'errpreset=conservative' not in line and circuit_initialization_parameters['simulation']['standard_parameters']['conservative']=='YES':
				line=line[:-1]+' errpreset=conservative \n'
		
		if write_check==1:
			s=s+line

	f.truncate(0)
	f.write(s)
	f.close()

#==========================================================================================================================


"""
===========================================================================================================================
------------------------------------ SPECTRE RUNNING FUNCTIONS ------------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will perform simulation for Basic Parameters
def write_extract_basic(circuit_initialization_parameters):
	
	# Writing the tcsh file for Basic Analysis
	sp.write_tcsh_file(circuit_initialization_parameters,'basic')

	# Writing the simulation parameters
	write_simulation_parameters(circuit_initialization_parameters,'basic')

	# Running netlist file
	sp.run_file(circuit_initialization_parameters)

	# Extracting the Basic Parameters
	basic_extracted_parameters=extract_basic_parameters(circuit_initialization_parameters)
	
	return basic_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will perform simulation for IIP3 Parameters
def write_extract_iip3(circuit_initialization_parameters):
	
	if circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']=='basic':
		
		# Writing the tcsh file for Basic Analysis
		sp.write_tcsh_file(circuit_initialization_parameters,'iip3')
		
		pin=circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']
		circuit_initialization_parameters['simulation']['netlist_parameters']['pin']=pin
		
		# Writing the simulation parameters
		write_simulation_parameters(circuit_initialization_parameters,'iip3')

		# Running netlist file
		sp.run_file(circuit_initialization_parameters)

		# Extracting Vout Magnitude
		file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.raw/hb_test.fd.qpss_hb'
		vout_fund_mag,vout_im3_mag=extract_vout_magnitude(file_name,circuit_initialization_parameters)

		# Calculating the iip3
		iip3_extracted_parameters={}
		iip3_extracted_parameters['iip3_dbm'],iip3_extracted_parameters['iip3_fund'],iip3_extracted_parameters['iip3_im3'],iip3_extracted_parameters['iip3_pin']=sp.calculate_iip3_single_point(vout_fund_mag,vout_im3_mag,pin)

	else:

		pin_start=circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']
		pin_stop=circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']
		pin_points=circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']

		pin=np.linspace(pin_start,pin_stop,pin_points)
		
		vout_fund_mag=np.zeros(pin_points,dtype=float)
		vout_im3_mag=np.zeros(pin_points,dtype=float)

		for i in range(pin_points):
		
			circuit_initialization_parameters['simulation']['netlist_parameters']['pin']=pin[i]
			
			# Writing the simulation parameters
			write_simulation_parameters(circuit_initialization_parameters)

			# Writing the tcsh file for Basic Analysis
			sp.write_tcsh_file(circuit_initialization_parameters,'iip3')

			# Running netlist file
			sp.run_file(circuit_initialization_parameters)

			# Extracting Vout Magnitude
			file_name=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.raw/hb_test.fd.qpss_hb'
			vout_fund_mag[i],vout_im3_mag[i]=extract_vout_magnitude(file_name,circuit_initialization_parameters)

		iip3_extracted_parameters={}
		n_pin=circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']
		n_points=circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']
		iip3_extracted_parameters['iip3_dbm'],iip3_extracted_parameters['iip3_im3_intercept'],iip3_extracted_parameters['iip3_im3_slope'],iip3_extracted_parameters['iip3_fund_intercept'],
		iip3_extracted_parameters['iip3_fund_slope']=sp.calculate_iip3_multiple_points(n_pin,n_points,vout_fund_mag,vout_im3_mag,pin)

	return iip3_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
def write_extract_single(i,circuit_parameters,circuit_initialization_parameters):
	
	# Writing to netlist file
	write_circuit_parameters(circuit_parameters,circuit_initialization_parameters)

	# Extracting the Basic Parameters
	basic_extracted_parameters=write_extract_basic(circuit_initialization_parameters)

	# Extracting the IIP3 Parameters
	iip3_extracted_parameters=write_extract_iip3(circuit_initialization_parameters)

	# Extracting Parameters from output files
	extracted_parameters=basic_extracted_parameters.copy()
	for param_name in iip3_extracted_parameters:
		extracted_parameters[param_name]=iip3_extracted_parameters[param_name]
	
	print(i,extracted_parameters)

	return (i,extracted_parameters)

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters for a single process
def write_extract(circuit_parameters,circuit_initialization_parameters):
	
	pool=mp.Pool(7)

	# Getting the values of frequency and range
	f_list=circuit_initialization_parameters['simulation']['standard_parameters']['f_list']
	n_freq=len(f_list)

	# Getting the different processes
	process_list=circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']
	n_process=len(process_list)

	# Getting the temperature list
	temp_list=circuit_initialization_parameters['simulation']['standard_parameters']['temp_list']
	n_temp=len(temp_list)

	# Getting the total number of runs
	n_runs=n_freq*n_process*n_temp

	# Creating new circuit parameter files
	circuit_parameters_run={}
	for i in range(n_runs):
		circuit_parameters_run[i]=circuit_parameters.copy()
		
	# Creating new circuit initialization parameters
	circuit_initialization_parameters_run={}
	for i in range(n_runs):
		i_freq,i_process,i_temp=sp.get_iteration(i,n_freq,n_process,n_temp)
		circuit_initialization_parameters_run[i]={}
		circuit_initialization_parameters_run[i]=copy.deepcopy(circuit_initialization_parameters)
		
		# Creating netlist directory
		netlist_folder=circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['directory']
		netlist_path=netlist_folder+'T'+str(i)+'/'
		if not os.path.exists(netlist_path):
			shutil.copytree(netlist_folder+'T_extra/',netlist_path)

		# Creating spectre run directory
		spectre_folder=circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['tcsh']+'Spectre_Run/'
		spectre_path=spectre_folder+'T'+str(i)+'/'
		if not os.path.exists(spectre_path):
			shutil.copytree(spectre_folder+'T_extra/',spectre_path)

		circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['directory']=netlist_path
		circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['tcsh']=spectre_path+'spectre_run.tcsh'
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['fund_1']=f_list[i_freq]
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['fund_2']=f_list[i_freq]+1e6
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['process_corner']=process_list[i_process]
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['cir_temp']=temp_list[i_temp]

	# Creating processes
	results_async=[pool.apply_async(write_extract_single,args=(i,circuit_parameters_run[i],circuit_initialization_parameters_run[i])) for i in range(n_runs)]

	extracted_parameters_combined={}
	for r in results_async:
		(i,extracted_parameters)=r.get()
		extracted_parameters_combined[i]=extracted_parameters
	
	extracted_parameters_split=sp.split_extracted_parameters(extracted_parameters_combined,f_list,process_list,temp_list)
	
	final_extracted_parameters=get_final_extracted_parameters(extracted_parameters_split,f_list,process_list,temp_list)
	
	pool.close()
	pool.join()

	return final_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters for a single process
def write_extract_multiple_circuits(circuit_parameters_dict,circuit_initialization_parameters):
	
	pool=mp.Pool(8)

	# Getting the values of frequency and range
	f_list=circuit_initialization_parameters['simulation']['standard_parameters']['f_list']
	n_freq=len(f_list)

	# Getting the different processes
	process_list=circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']
	n_process=len(process_list)

	# Getting the temperature list
	temp_list=circuit_initialization_parameters['simulation']['standard_parameters']['temp_list']
	n_temp=len(temp_list)

	# Getting the circuit parameters list
	n_circuits=len(circuit_parameters_dict)

	# Getting the total number of runs
	n_runs=n_circuits*n_freq*n_process*n_temp

	# Creating new circuit parameter files
	circuit_parameters_run={}
	for i in range(n_runs):
		circuit_parameters_run[i]=circuit_parameters_dict[i%n_circuits].copy()

		print('~~~~~~~~~~~~~~~~~~~~~~~~~ CIRCUIT PARAMETER RUN '+str(i)+'~~~~~~~~~~~~~~~~')
		
	# Creating new circuit initialization parameters
	circuit_initialization_parameters_run={}
	for i in range(n_runs):
		i_freq,i_process,i_temp=sp.get_iteration(i//n_circuits,n_freq,n_process,n_temp)
		circuit_initialization_parameters_run[i]={}
		circuit_initialization_parameters_run[i]=copy.deepcopy(circuit_initialization_parameters)
		
		# Creating netlist directory
		netlist_folder=circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['directory']
		netlist_path=netlist_folder+'T'+str(i)+'/'
		if not os.path.exists(netlist_path):
			shutil.copytree(netlist_folder+'T_extra/',netlist_path)

		# Creating spectre run directory
		spectre_folder=circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['tcsh']+'Spectre_Run/'
		spectre_path=spectre_folder+'T'+str(i)+'/'
		if not os.path.exists(spectre_path):
			shutil.copytree(spectre_folder+'T_extra/',spectre_path)

		circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['directory']=netlist_path
		circuit_initialization_parameters_run[i]['simulation']['standard_parameters']['tcsh']=spectre_path+'spectre_run.tcsh'
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['fund_1']=f_list[i_freq]
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['fund_2']=f_list[i_freq]+1e6
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['process_corner']=process_list[i_process]
		circuit_initialization_parameters_run[i]['simulation']['netlist_parameters']['cir_temp']=temp_list[i_temp]
		
	# Creating processes
	results_async=[pool.apply_async(write_extract_single,args=(i,circuit_parameters_run[i],circuit_initialization_parameters_run[i])) for i in range(n_runs)]

	extracted_parameters_combined={}
	for r in results_async:
		(i,extracted_parameters)=r.get()
		extracted_parameters_combined[i]=extracted_parameters
		
	extracted_parameters_split=sp.split_extracted_parameters_multiple(extracted_parameters_combined,f_list,process_list,temp_list,n_circuits)
	
	final_extracted_parameters={}
	for i in range(n_circuits):
		final_extracted_parameters[i]=get_final_extracted_parameters(extracted_parameters_split[i],f_list,process_list,temp_list)
	
	pool.close()
	pool.join()

	return final_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters
def get_final_extracted_parameters(extracted_parameters_split,f_list,process_list,temp_list):
	
	print('\n\n\nExtracted Parameters Split\n')
	print(extracted_parameters_split)
	extracted_parameters_frequency=get_final_extracted_parameters_frequency(extracted_parameters_split,f_list,process_list,temp_list)
	
	print('\n\n\nExtracted Parameters Frequency\n')
	print(extracted_parameters_frequency)
	extracted_parameters_process=get_final_extracted_parameters_process(extracted_parameters_frequency,process_list,temp_list)
	
	print('\n\n\nExtracted Parameters Process\n')
	print(extracted_parameters_process)
	final_extracted_parameters=get_final_extracted_parameters_temperature(extracted_parameters_process,temp_list)
	
	print('\n\n\nFinal Extracted Parameters Split\n')
	print(final_extracted_parameters)
	return final_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run spectre and extract the output parameters
def get_final_extracted_parameters_frequency(extracted_parameters_split,f_list,process_list,temp_list):
	
	if len(f_list)==1:
		extracted_parameters_frequency={}
		for temp in temp_list:
			extracted_parameters_frequency[temp]={}
			for process in process_list:
				extracted_parameters_frequency[temp][process]={}
				for freq in extracted_parameters_split[temp][process]:
					for param in extracted_parameters_split[temp][process][freq]:
						extracted_parameters_frequency[temp][process][param]=extracted_parameters_split[temp][process][freq][param]
	
	else:
		# Getting the least, middle, and maximum frequency
		f_len=len(f_list)
		small_frequency=f_list[0]
		mid_frequency=f_list[(f_len-1)//2]
		large_frequency=f_list[f_len-1]

		# Selecting which parameters are DC, and which parameter to select for AC among different frequencies
		extracted_parameters_select={
			'v_source':'dc',
			'i_source':'dc',
			'p_source':'dc',
			
			'vg':'dc',
			'vd':'dc',
			'vs':'dc',

			'Io':'dc',
			'gm1':'dc',
			'gds1':'dc',
			'vt':'dc',
			'vdsat':'dc',
			'cgs1':'dc',
			'cgd1':'dc',
			
			'freq':'mid',
			's12_db':'max',
			's21_db':'max',
			's22_db':'max',
			'k':'min',
			'nf_db':'max',
			'iip3_dbm':'min'
		}

		# First, we will iterate among different temperature and process lists
		extracted_parameters_frequency={}
		for temp in temp_list:
			extracted_parameters_frequency[temp]={}
			for process in process_list:
				extracted_parameters_frequency[temp][process]={}

				# Now, we have a given process and temperature ; We need to find the extracted parameters for this set
				final_extracted_parameters={}
				extracted_parameters_combined=extracted_parameters_split[temp][process]

				# Getting the values for all frequencies in case it is AC
				for param in extracted_parameters_combined[mid_frequency]:
					if param in extracted_parameters_select:
						if extracted_parameters_select[param]=='dc':
							continue
					for freq in f_list:
						final_extracted_parameters[str(freq)+'_'+param]=extracted_parameters_combined[freq][param]
				
				# Getting the min or mid or max parameter values for the best values
				for param in extracted_parameters_select:

					# Case I : DC Parameter or AC Parameter with mid value
					if extracted_parameters_select[param]=='mid' or extracted_parameters_select[param]=='dc':
						final_extracted_parameters[param]=extracted_parameters_combined[mid_frequency][param]

					# Case II : AC Parameter with minimum value
					elif extracted_parameters_select[param]=='min':
						param_array=[]
						for i in extracted_parameters_combined:
							param_array.append(extracted_parameters_combined[i][param])
						final_extracted_parameters[param]=min(param_array)
					
					# Case III : AC Parameter with maximum value
					else:
						param_array=[]
						for i in extracted_parameters_combined:
							param_array.append(extracted_parameters_combined[i][param])
						final_extracted_parameters[param]=max(param_array)
				
				# Calculating the value of gain
				gain_array=[]
				gain_phase_array=[]
				for i in extracted_parameters_combined:
					gain_array.append(extracted_parameters_combined[i]['gain_db'])
					gain_phase_array.append(extracted_parameters_combined[i]['gain_phase'])
				gain_min=min(gain_array)
				gain_index=f_list[gain_array.index(gain_min)]
				final_extracted_parameters['gain_db']=gain_min
				final_extracted_parameters['gain_phase']=extracted_parameters_combined[gain_index]['gain_phase']
				
				# Calculating the value of s11
				s11_array=[]
				ZR_array=[]
				ZI_array=[]
				for i in extracted_parameters_combined:
					s11_array.append(extracted_parameters_combined[i]['s11_db'])
					ZR_array.append(extracted_parameters_combined[i]['Zin_R'])
					ZI_array.append(extracted_parameters_combined[i]['Zin_I'])
				s11_max=max(s11_array)
				s11_index=f_list[s11_array.index(s11_max)]
				final_extracted_parameters['s11_db']=s11_max
				final_extracted_parameters['Zin_R']=extracted_parameters_combined[s11_index]['Zin_R']
				final_extracted_parameters['Zin_I']=extracted_parameters_combined[s11_index]['Zin_I']
				
				# Getting the iip3 values
				iip3_array_list=['iip3_im3_intercept','iip3_im3_slope','iip3_fund_intercept','iip3_fund_slope','iip3_fund','iip3_im3','iip3_pin']
				for param in iip3_array_list:
					if param in extracted_parameters_combined[mid_frequency]:
						final_extracted_parameters[param]=extracted_parameters_combined[mid_frequency][param]

				extracted_parameters_frequency[temp][process]=final_extracted_parameters.copy()

	return extracted_parameters_frequency

#-----------------------------------------------------------------------------------------------
# This function will combine the extracted_parameters
def get_final_extracted_parameters_process(extracted_parameters_frequency,process_list,temp_list):
	
	# Getting the centre process
	process_len=len(process_list)
	middle_process=process_list[(process_len-1)//2]

	# Getting the values for all process for all temperature
	extracted_parameters_process={}
	for temp in temp_list:
		extracted_parameters_process[temp]={}

		if len(process_list)==1:
			for process in process_list:
				for param_name in extracted_parameters_frequency[temp][process]:
					extracted_parameters_process[temp][param_name]=extracted_parameters_frequency[temp][process][param_name]

		else:
			for process in process_list:
				for param_name in extracted_parameters_frequency[temp][process]:
					extracted_parameters_process[temp][process+'_'+param_name]=extracted_parameters_frequency[temp][process][param_name]
	
	if len(process_list)!=1:
		extracted_parameters_select={
			'v_source':'mid',
			'i_source':'mid',
			'p_source':'mid',
			
			'vg':'mid',
			'vd':'mid',
			'vs':'mid',

			'Io':'mid',
			'gm1':'mid',
			'gds1':'mid',
			'vt':'mid',
			'vdsat':'mid',
			'cgs1':'mid',
			'cgd1':'mid',
			
			'freq':'mid',
			's11_db':'max',
			's12_db':'mid',
			's21_db':'mid',
			's22_db':'mid',
			'k':'mid',
			'nf_db':'max',
			'iip3_dbm':'min',

			'gain_db':'min',
			'Zin_R':'mid',
			'Zin_I':'mid'
		}

		# Choosing the best process among the different temperatures
		for temp in temp_list:
			
			# Getting the min or mid or max parameter values for the best values
			for param in extracted_parameters_select:

				# Case I : Minimum value
				if extracted_parameters_select[param]=='min':
					param_array=[]
					for process in process_list:
						param_array.append(extracted_parameters_frequency[temp][process][param])
					extracted_parameters_process[temp][param]=min(param_array)
				
				# Case II : Maximum value
				elif extracted_parameters_select[param]=='max':
					param_array=[]
					for process in process_list:
						param_array.append(extracted_parameters_frequency[temp][process][param])
					extracted_parameters_process[temp][param]=max(param_array)
				
				# Case III : Middle value ( for typical corner )
				else:
					extracted_parameters_process[temp][param]=extracted_parameters_frequency[temp][middle_process][param]
	
	return extracted_parameters_process

#-----------------------------------------------------------------------------------------------
# This function will combine the extracted_parameters
def get_final_extracted_parameters_temperature(extracted_parameters_process,temp_list):
	
	# Getting the centre temperature
	temp_len=len(temp_list)
	middle_temp=temp_list[(temp_len-1)//2]

	# Getting the values for all temperatures
	extracted_parameters={}

	if len(temp_list)==1:
		for temp in extracted_parameters_process:
			for param_name in extracted_parameters_process[temp]:
				extracted_parameters[param_name]=extracted_parameters_process[temp][param_name]

	else:
		for temp in extracted_parameters_process:
			for param_name in extracted_parameters_process[temp]:
				extracted_parameters[str(temp)+'_'+param_name]=extracted_parameters_process[temp][param_name]
	
	if len(temp_list)!=1:
		extracted_parameters_select={
			'v_source':'mid',
			'i_source':'mid',
			'p_source':'mid',
			
			'vg':'mid',
			'vd':'mid',
			'vs':'mid',

			'Io':'mid',
			'gm1':'mid',
			'gds1':'mid',
			'vt':'mid',
			'vdsat':'mid',
			'cgs1':'mid',
			'cgd1':'mid',
			
			'freq':'mid',
			's11_db':'max',
			's12_db':'mid',
			's21_db':'mid',
			's22_db':'mid',
			'k':'mid',
			'nf_db':'max',
			'iip3_dbm':'min',

			'gain_db':'min',
			'Zin_R':'mid',
			'Zin_I':'mid'
		}

		# Getting the min or mid or max parameter values for the best values
		for param in extracted_parameters_select:
			
			# Case I : Minimum value
			if extracted_parameters_select[param]=='min':
				param_array=[]
				for temp in temp_list:
					param_array.append(extracted_parameters_process[temp][param])
				extracted_parameters[param]=min(param_array)
			
			# Case II : Maximum value
			elif extracted_parameters_select[param]=='max':
				param_array=[]
				for temp in temp_list:
					param_array.append(extracted_parameters_process[temp][param])
				extracted_parameters[param]=max(param_array)
			
			# Case III : Middle value ( for central temperature )
			else:
				extracted_parameters[param]=extracted_parameters_process[middle_temp][param]
		
	return extracted_parameters



#==========================================================================================================================


#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
def write_extract_old(circuit_parameters,circuit_initialization_parameters):

	# Getting the operating frequency
	f_operating=circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']
	circuit_initialization_parameters['simulation']['netlist_parameters']['fund_1']=f_operating
	circuit_initialization_parameters['simulation']['netlist_parameters']['fund_2']=f_operating+1e6
	
	# Writing to netlist file
	write_circuit_parameters(circuit_parameters,circuit_initialization_parameters)

	# Extracting the Basic Parameters
	basic_extracted_parameters=write_extract_basic(circuit_initialization_parameters)

	# Extracting the IIP3 Parameters
	iip3_extracted_parameters=write_extract_iip3(circuit_initialization_parameters)

	# Extracting Parameters from output files
	extracted_parameters=basic_extracted_parameters.copy()
	for param_name in iip3_extracted_parameters:
		extracted_parameters[param_name]=iip3_extracted_parameters[param_name]
	
	return extracted_parameters

#==========================================================================================================================
