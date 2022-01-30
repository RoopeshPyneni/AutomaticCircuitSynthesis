#==========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will contain the functions to write, run, and read from the spectre files for CS LNA
"""

#==========================================================================================================================
import numpy as np
import fileinput
import multiprocessing as mp
import CS_LNA.pre_optimization as pr # type: ignore
import spectre_common as sp # type: ignore
import copy
import sys

"""
===========================================================================================================================
------------------------------------- CIRCUIT CLASS -----------------------------------------------------------------------
"""

#--------------------------------------------------------------------------------------------------------------------------
# Creating a class for the circuit 
class Circuit():
	def __init__(self,circuit_initialization_parameters):
		self.circuit_parameters={}
		self.extracted_parameters={}
		self.simulation_parameters={}
		self.circuit_initialization_parameters=circuit_initialization_parameters
		self.mos_parameters=sp.calculate_mos_parameters(self.circuit_initialization_parameters)

		# Getting the circuit name
		if self.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_resistor':
			self.circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_r'
			self.circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_r'
		
		elif self.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='series':
			self.circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_series'
			self.circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_series'

		elif self.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_capacitor':
			self.circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters_rc'
			self.circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb_rc'

		elif self.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='ideal':
			self.circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters'
			self.circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb'
		
		else:
			sys.exit()
	
	def run_circuit(self):
		self.extracted_parameters=write_extract(self.circuit_parameters,self.circuit_initialization_parameters)

	def update_circuit(self,circuit_parameters):
		self.circuit_parameters=circuit_parameters
		if self.circuit_parameters['Rk']>=1.00:
			self.circuit_parameters['Rk']=0.90
		self.extracted_parameters=write_extract(circuit_parameters,self.circuit_initialization_parameters)
	
	def update_circuit_parameters(self,circuit_parameters):
		self.circuit_parameters=circuit_parameters
		if self.circuit_parameters['Rk']>=1.00:
			self.circuit_parameters['Rk']=0.90

	def update_simulation_parameters(self,simulation_parameters):
		if 'netlist_parameters' in simulation_parameters:
			for param_name in simulation_parameters['netlist_parameters']:
				self.circuit_initialization_parameters['simulation']['netlist_parameters'][param_name]=simulation_parameters['netlist_parameters'][param_name]
		
		if 'standard_parameters' in simulation_parameters:
			for param_name in simulation_parameters['standard_parameters']:
				self.circuit_initialization_parameters['simulation']['standard_parameters'][param_name]=simulation_parameters['standard_parameters'][param_name]
	
	def write_simulation_parameters(self):
		write_simulation_parameters(self.circuit_initialization_parameters)
	
	def update_temp(self,temp):
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=temp
	
	def reset_temp(self):
		self.circuit_initialization_parameters['simulation']['netlist_parameters']['cir_temp']=self.circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']	

	def calculate_iip3(self,n_pin,n_points,vout_fund_mag,vout_im3_mag,pin):
		return sp.calculate_iip3_multiple_points(n_pin,n_points,vout_fund_mag,vout_im3_mag,pin)

	
	#-----------------------------------------------------------------------------------------------
	# Optimization Functions

	#-----------------------------------------------------------------------------------------------
	# This function calculates the loss for Io Optimization
	def calc_loss(self,output_conditions,loss_weights):
		
		# Extracted Values
		gain=self.extracted_parameters['gain_db']
		iip3=self.extracted_parameters['iip3_dbm']
		s11=self.extracted_parameters['s11_db']
		nf=self.extracted_parameters['nf_db']
		Io=self.extracted_parameters['Io']

		gain_0=self.extracted_parameters['0_gain_db']
		gain_2=self.extracted_parameters['2_gain_db']
		
		# Reference Values
		gain_ref=output_conditions['gain_db']
		iip3_ref=output_conditions['iip3_dbm']
		s11_ref=output_conditions['s11_db']
		nf_ref=output_conditions['nf_db']

		gain_delta_ref=output_conditions['gain_delta']
		
		#Defining the weights to calculate Loss
		A1=loss_weights['gain_db']	# Weight for gain
		A2=loss_weights['iip3_dbm']	# Weight for iip3
		A3=loss_weights['s11_db']	# Weight for s11
		A4=loss_weights['nf_db']	# Weight for nf
		A5=loss_weights['Io']	# Weight for Io
		
		A6=loss_weights['gain_delta']
		
		# Calculating Loss
		loss_gain=A1*sp.ramp_func(gain_ref-gain)
		loss_iip3=A2*sp.ramp_func(iip3_ref-iip3)
		loss_s11=A3*sp.ramp_func(s11-s11_ref)
		loss_nf=A4*sp.ramp_func(nf-nf_ref)
		loss_Io=A5*Io
		loss_gain_delta=A6*sp.ramp_func(gain_0-gain)+A6*sp.ramp_func(gain-gain_delta_ref-gain_0)+A6*sp.ramp_func(gain_2-gain)+A6*sp.ramp_func(gain-gain_delta_ref-gain_2)
		loss=loss_gain+loss_iip3+loss_s11+loss_nf+loss_Io+loss_gain_delta
		loss_dict={'loss':loss,'loss_gain':loss_gain,'loss_iip3':loss_iip3,'loss_s11':loss_s11,'loss_nf':loss_nf,'loss_Io':loss_Io,'loss_gain_delta':loss_gain_delta}
		
		return loss_dict
	
	#-----------------------------------------------------------------------------------------------
	# This function updates the values of circuit parameters by trying to minimize loss
	def update_circuit_parameters(self,circuit_parameters_slope,check_loss,optimization_input_parameters,run_number):
		
		alpha=optimization_input_parameters['optimization'][run_number]['alpha']['value']

		# Calculating the value to update each parameter with
		for param_name in circuit_parameters_slope:
			
			# Calculating the Increment Value
			if check_loss==-1:
				change=circuit_parameters_slope[param_name]['loss']*(self.circuit_parameters[param_name]**2)*alpha
			elif check_loss==1:
				change=circuit_parameters_slope[param_name]['loss_Io']*(self.circuit_parameters[param_name]**2)*alpha
			else:
				change=(circuit_parameters_slope[param_name]['loss']-circuit_parameters_slope[param_name]['loss_Io'])
				change=change*(self.circuit_parameters[param_name]**2)*alpha
		
		
			# Checking if the parameter is updated by a large value
			change_limit=0.25 # If the incremented value is more than +- change_limit*parameter_name, then we will limit the change
			if change>change_limit*self.circuit_parameters[param_name]:
				change=change_limit*self.circuit_parameters[param_name]
			if change<-1*change_limit*self.circuit_parameters[param_name]:
				change=-1*change_limit*self.circuit_parameters[param_name]
			
			# Updating circuit_parameters
			self.circuit_parameters[param_name]=self.circuit_parameters[param_name]-change
			
	#-----------------------------------------------------------------------------------------------
	# This function will check the loss of gain, iip3, nf, and s11
	def calc_check_loss(self,loss_iter,i,loss_type):

		if loss_type==0:
			if loss_iter[i-1]['loss']==loss_iter[i-1]['loss_Io']:
				check_loss=1
			else:
				check_loss=0
					
		elif loss_type==1:
			check_loss=-1
			
		return check_loss

	#---------------------------------------------------------------------------------------------------------------------------
	# Function to check the best solution
	def check_best_solution(self,optimization_results,loss_max):

		# Defining some values
		n_iter=optimization_results['n_iter']
		iter_min=0
		
		zero_loss_array=['loss_s11','loss_gain','loss_iip3','loss_nf','loss_gain_delta']
		minimize_loss_array=['loss_Io']

		loss_Io_min=sum([optimization_results['loss_iter'][0][key] for key in minimize_loss_array])

		if sum([optimization_results['loss_iter'][0][key] for key in zero_loss_array])>loss_max:
			flag=0
		else:
			flag=1
		
		for i in range(1,n_iter):
			if sum([optimization_results['loss_iter'][i][key] for key in zero_loss_array])>loss_max:
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

	#---------------------------------------------------------------------------------------------------------------------------
	# Function to perform pre optimization
	def pre_optimization(self,optimization_input_parameters,timing_results):
		pr.pre_optimization(self,optimization_input_parameters,timing_results)


"""
===========================================================================================================================
------------------------------------- EXTRACTION FUNCTION -----------------------------------------------------------------
"""

#==========================================================================================================================
#------------------------------------ Basic File Extraction Functions -----------------------------------------------------

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
	extracted_parameters['v_source']=np.absolute(sp.valueE_to_value(lines[1]))
	extracted_parameters['i_source']=np.absolute(sp.valueE_to_value(lines[2]))
	extracted_parameters['p_source']=extracted_parameters['i_source']*extracted_parameters['v_source']

	extracted_parameters['vg1']=sp.valueE_to_value(lines[3])
	extracted_parameters['vd1']=sp.valueE_to_value(lines[4])
	extracted_parameters['vg2']=sp.valueE_to_value(lines[5])
	
	extracted_parameters['Io']=sp.valueE_to_value(lines[6])
	extracted_parameters['gm1']=sp.valueE_to_value(lines[7])
	extracted_parameters['gds1']=sp.valueE_to_value(lines[8])
	extracted_parameters['vth1']=sp.valueE_to_value(lines[9])
	extracted_parameters['vds1']=sp.valueE_to_value(lines[10])
	extracted_parameters['vdsat1']=sp.valueE_to_value(lines[11])
	extracted_parameters['cgs1']=np.absolute(sp.valueE_to_value(lines[12]))
	extracted_parameters['cgd1']=np.absolute(sp.valueE_to_value(lines[13]))
	extracted_parameters['region1']=sp.valueE_to_value(lines[14])
	if extracted_parameters['vds1']-extracted_parameters['vdsat1']>0.1:
		extracted_parameters['check_vd1']=1
	else:
		extracted_parameters['check_vd1']=0

	extracted_parameters['vth2']=sp.valueE_to_value(lines[15])
	extracted_parameters['vds2']=sp.valueE_to_value(lines[16])
	extracted_parameters['vdsat2']=sp.valueE_to_value(lines[17])
	extracted_parameters['cgs2']=np.absolute(sp.valueE_to_value(lines[18]))
	extracted_parameters['cgd2']=np.absolute(sp.valueE_to_value(lines[19]))
	extracted_parameters['region2']=sp.valueE_to_value(lines[20])
	if extracted_parameters['vds2']-extracted_parameters['vdsat2']>0.1:
		extracted_parameters['check_vd2']=1
	else:
		extracted_parameters['check_vd2']=0

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
	extracted_parameters['vgs_ac']=sp.valueE_to_value(lines[5])
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
	#lines=lines[12:]
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
#------------------------------------ IIP3 Extraction Functions -----------------------------------------------------------

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
------------------------------------- FILE WRITE FUNCTIONS ----------------------------------------------------------------
"""

#-----------------------------------------------------------------      
# Function that converts input parameter dictionary to writing dictionary
def dict_convert(circuit_parameters,circuit_initialization_parameters):
	write_dict={}
	
	# param_names in write_dict will contain the name of the parameters as it is written in the .scs file
	cir_writing_dict={	
		'wid':'W',
		'cur0':'Io',
		'res_b':'Rb',
		'ind_d':'Ld',
		'ind_g':'Lg',
		'ind_s':'Ls',
		'cap_s':'Cs',
		'cap_g':'Cg',
		'cap_d':'Cd',
		'res_sum':'Rsum',
		'res_k':'Rk',
		'res_l':'Rl'
	}

	# Checking for wrong values
	if circuit_parameters['Rk']>=1.0:
		circuit_parameters['Rk']=0.95
	if circuit_parameters['Rk']<=0.0:
		circuit_parameters['Rk']=0.05

	for param_name in cir_writing_dict:
		write_dict[param_name]=circuit_parameters[cir_writing_dict[param_name]]
	
	# Getting the value of resistances
	write_dict['res_g']=circuit_parameters['Lg']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi*50
	write_dict['res_d']=circuit_parameters['Ld']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi*15
	write_dict['res_ls']=circuit_parameters['Ls']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi*15

	if circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='series':
		write_dict['res_g']=circuit_parameters['Lg']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi/50
		write_dict['res_ls']=circuit_parameters['Ls']*circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']*2*np.pi/15

	# Calculating the number of fingers
	n_finger=int(circuit_parameters['W']/circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max'])+1
	write_dict['n_finger']=n_finger

	# Getting the width and length for TSMC Resistors
	write_dict['res_b_len'],write_dict['res_b_wid']=get_TSMC_resistor(circuit_parameters['Rb'])
	R1=circuit_parameters['Rsum']*(1-circuit_parameters['Rk'])
	R2=circuit_parameters['Rsum']*circuit_parameters['Rk']
	write_dict['res_1_len'],write_dict['res_1_wid']=get_TSMC_resistor(R1)
	write_dict['res_2_len'],write_dict['res_2_wid']=get_TSMC_resistor(R2)
	write_dict['res_l_len'],write_dict['res_l_wid']=get_TSMC_resistor(circuit_parameters['Rl'])

	# Getting the width, length, mf for Capacitors
	write_dict['wid_cap_g'],write_dict['len_cap_g']=calculate_MOS_capacitor(circuit_parameters['Cg'])
	write_dict['mf_cap_s']=1+int(circuit_parameters['Cs']*1e11)
	write_dict['mf_cap_d']=1+int(circuit_parameters['Cd']/5.63e-15)

	return write_dict

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
# Function that modifies the .scs file
def write_circuit_parameters(circuit_parameters,circuit_initialization_parameters):
	
	# We will convert the circuit parameters to write_dict
	write_dict=dict_convert(circuit_parameters,circuit_initialization_parameters)
	
	# Getting the filenames
	filename1=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
	filename2=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.scs'

	# We will write the new values to the Basic Circuit
	f=open(filename1,'r+')
	s=''
	for line in fileinput.input(filename1):
		for param_name in write_dict:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,sp.print_param(param_name,write_dict[param_name]))	# Replacing the parameter in the .scs file
		s=s+line
	f.truncate(0)
	f.write(s)
	f.close()

	# We will write the new values to the IIP3 Circuit
	f=open(filename2,'r+')
	s=''
	for line in fileinput.input(filename2):
		for param_name in write_dict:
			if "parameters "+param_name+'=' in line:	# Checking for a particular parameter in the .scs file
				line=line.replace(line,sp.print_param(param_name,write_dict[param_name]))	# Replacing the parameter in the .scs file
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
def write_simulation_parameters(circuit_initialization_parameters):
	
	# Adding simulation_parameters to write_dict
	write_dict={}
	for param_name in circuit_initialization_parameters['simulation']['netlist_parameters']:
		write_dict[param_name]=circuit_initialization_parameters['simulation']['netlist_parameters'][param_name]
	process_corner=circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']
	
	write_dict['len']=circuit_initialization_parameters['MOS']['Lmin']
	write_dict['v_dd']=circuit_initialization_parameters['MOS']['Vdd']

	# Getting the filenames
	filename1=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']+'/circ.scs'
	filename2=circuit_initialization_parameters['simulation']['standard_parameters']['directory']+circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']+'/circ.scs'

	# Writing the simulation parameters to Basic File
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

	# Writing the simulation parameters to IIP3 File
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

#==========================================================================================================================


"""
===========================================================================================================================
------------------------------------- SPECTRE RUNNING FUNCTIONS -----------------------------------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will perform simulation for Basic Parameters
def write_extract_basic(circuit_initialization_parameters):
	
	# Writing the tcsh file for Basic Analysis
	sp.write_tcsh_file(circuit_initialization_parameters,'basic')

	# Writing the simulation parameters
	write_simulation_parameters(circuit_initialization_parameters)

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
		write_simulation_parameters(circuit_initialization_parameters)

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

	return (i,extracted_parameters)

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
def get_final_extracted_parameters(extracted_parameters_combined):
	
	final_extracted_parameters={}

	extracted_parameters_select={
		'v_source':'dc',
		'i_source':'dc',
		'p_source':'dc',
		
		'vg1':'dc',
		'vd1':'dc',
		'vg2':'dc',

		'Io':'dc',
		'gm1':'dc',
		'gds1':'dc',
		'vth1':'dc',
		'vds1':'dc',
		'vdsat1':'dc',
		'cgs1':'dc',
		'cgd1':'dc',
		'region1':'dc',
		'check_vd1':'dc',

		'vth2':'midcd',
		'vds2':'dc',
		'vdsat2':'dc',
		'cgs2':'dc',
		'cgd2':'dc',
		'region2':'dc',
		'check_vd2':'dc',

		'freq':'mid',
		'vgs_ac':'mid',
		's12_db':'max',
		's21_db':'max',
		's22_db':'max',
		'k':'min',
		'nf_db':'max',
		'iip3_dbm':'min'
	}

	# Getting the values for all three frequencies
	for param in extracted_parameters_combined[0]:
		if param in extracted_parameters_select:
			if extracted_parameters_select[param]=='dc':
				continue
		for i in range(3):
			final_extracted_parameters[str(i)+'_'+param]=extracted_parameters_combined[i][param]
	
	# Getting the min or mid or max parameter values for the best values
	for param in extracted_parameters_select:
		if extracted_parameters_select[param]=='mid' or extracted_parameters_select[param]=='dc':
			final_extracted_parameters[param]=extracted_parameters_combined[1][param]
		elif extracted_parameters_select[param]=='min':
			param_array=[]
			for i in extracted_parameters_combined:
				param_array.append(extracted_parameters_combined[i][param])
			final_extracted_parameters[param]=min(param_array)
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
	gain_index=gain_array.index(gain_min)
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
	s11_index=s11_array.index(s11_max)
	final_extracted_parameters['s11_db']=s11_max
	final_extracted_parameters['Zin_R']=extracted_parameters_combined[s11_index]['Zin_R']
	final_extracted_parameters['Zin_I']=extracted_parameters_combined[s11_index]['Zin_I']

	# Getting the iip3 values
	iip3_array_list=['iip3_im3_intercept','iip3_im3_slope','iip3_fund_intercept','iip3_fund_slope','iip3_fund','iip3_im3','iip3_pin']
	for param in iip3_array_list:
		if param in extracted_parameters_combined[1]:
			final_extracted_parameters[param]=extracted_parameters_combined[1][param]

	return final_extracted_parameters

#-----------------------------------------------------------------------------------------------
# This function will write the circuit parameters, run Eldo and extract the output parameters
def write_extract(circuit_parameters,circuit_initialization_parameters):
	
	pool=mp.Pool()

	# Creating new circuit parameter files
	circuit_parameters_run={}
	circuit_parameters_run[0]=circuit_parameters.copy()
	circuit_parameters_run[1]=circuit_parameters.copy()
	circuit_parameters_run[2]=circuit_parameters.copy()
	circuit_initialization_parameters_run={}
	

	# Getting the values of frequency and range
	f_operating=circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']
	f_range=circuit_initialization_parameters['simulation']['standard_parameters']['f_range']

	# Creating new circuit initialization parameters
	circuit_initialization_parameters_run[0]={}
	circuit_initialization_parameters_run[0]=copy.deepcopy(circuit_initialization_parameters)
	circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['directory']=circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['directory']+'T1/'
	circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[0]['simulation']['standard_parameters']['tcsh']+'Spectre_Run/T1/spectre_run.tcsh'
	circuit_initialization_parameters_run[0]['simulation']['netlist_parameters']['fund_1']=f_operating-f_range
	circuit_initialization_parameters_run[0]['simulation']['netlist_parameters']['fund_2']=f_operating-f_range+1e6


	circuit_initialization_parameters_run[1]={}
	circuit_initialization_parameters_run[1]=copy.deepcopy(circuit_initialization_parameters)
	circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['directory']=circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['directory']+'T2/'
	circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[1]['simulation']['standard_parameters']['tcsh']+'Spectre_Run/T2/spectre_run.tcsh'
	circuit_initialization_parameters_run[1]['simulation']['netlist_parameters']['fund_1']=f_operating
	circuit_initialization_parameters_run[1]['simulation']['netlist_parameters']['fund_2']=f_operating+1e6
	

	circuit_initialization_parameters_run[2]={}
	circuit_initialization_parameters_run[2]=copy.deepcopy(circuit_initialization_parameters)
	circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['directory']=circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['directory']+'T3/'
	circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['tcsh']=circuit_initialization_parameters_run[2]['simulation']['standard_parameters']['tcsh']+'Spectre_Run/T3/spectre_run.tcsh'
	circuit_initialization_parameters_run[2]['simulation']['netlist_parameters']['fund_1']=f_operating+f_range
	circuit_initialization_parameters_run[2]['simulation']['netlist_parameters']['fund_2']=f_operating+f_range+1e6
		

	# Creating processes
	results_async=[pool.apply_async(write_extract_single,args=(i,circuit_parameters_run[i],circuit_initialization_parameters_run[i])) for i in range(3)]

	extracted_parameters_combined={}
	for r in results_async:
		(i,extracted_parameters)=r.get()
		extracted_parameters_combined[i]=extracted_parameters
		
	final_extracted_parameters=get_final_extracted_parameters(extracted_parameters_combined)
	
	pool.close()
	pool.join()

	return final_extracted_parameters


#==========================================================================================================================
