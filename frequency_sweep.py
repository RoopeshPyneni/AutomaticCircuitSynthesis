#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will run spectre and extract the parameters for a range of frequencies
"""

#===========================================================================================================================
#import CG_LNA.spectre as sp
import CS_LNA.spectre as sp
import numpy as np
import os
from matplotlib import pylab
from pylab import *


"""
===========================================================================================================================
------------------------------------ Other Functions ----------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the MOSFET parameters to the circuit_initialization_parameters dictionary
def get_mos_parameters(circuit_initialization_parameters,process_name):
	
	circuit_initialization_parameters['MOS']={}
	circuit_initialization_parameters['MOS']['Process']=process_name
	circuit_initialization_parameters['MOS']['filename']={}
	
	f=open('/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/MOS_Files/'+process_name+'.txt')
	lines=f.readlines()
	f.close()

	# Extracting values from the MOS File
	for i in range(len(lines)):
		line=lines[i][:-1]
		if line=='Vdd':
			circuit_initialization_parameters['MOS']['Vdd']=float(lines[i+1][:-1])
		elif line=='Lmin':
			circuit_initialization_parameters['MOS']['Lmin']=float(lines[i+1][:-1])
		elif line=='u0':
			circuit_initialization_parameters['MOS']['un']=float(lines[i+1][:-1])*1e-4
		elif line=='tox':
			circuit_initialization_parameters['MOS']['tox']=float(lines[i+1][:-1])
		elif line=='vth0':
			circuit_initialization_parameters['MOS']['vt']=float(lines[i+1][:-1])
		elif line=='tt_file':
			circuit_initialization_parameters['MOS']['filename']['tt']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['tt']+=lines[j]
				j+=1
		elif line=='ff_file':
			circuit_initialization_parameters['MOS']['filename']['ff']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['ff']+=lines[j]
				j+=1
		elif line=='ss_file':
			circuit_initialization_parameters['MOS']['filename']['ss']=''
			j=i+1
			while lines[j][:-1]!='':
				circuit_initialization_parameters['MOS']['filename']['ss']+=lines[j]
				j+=1
                
	# Calculating Cox
	eo=8.85*1e-12
	er=3.9
	circuit_initialization_parameters['MOS']['cox']=eo*er/circuit_initialization_parameters['MOS']['tox']

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions_CG_LNA(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}
	circuit_initialization_parameters['simulation']['directory']='/home/ee18b028/cadence_project/lna1/'
	circuit_initialization_parameters['simulation']['basic_circuit']='basic_parameters_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['iip3_circuit']='iip3_hb_tsmc_65_rcm'
	circuit_initialization_parameters['simulation']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/spectre_run.tcsh'
	circuit_initialization_parameters['simulation']['iip3_type']='basic'		# 'basic' or 'advanced' 

	circuit_initialization_parameters['simulation']['std_temp']=27
	circuit_initialization_parameters['simulation']['pin_fixed']=-65
	circuit_initialization_parameters['simulation']['pin_start']=-70
	circuit_initialization_parameters['simulation']['pin_stop']=-40
	circuit_initialization_parameters['simulation']['pin_points']=6
	circuit_initialization_parameters['simulation']['iip3_calc_points']=3
	circuit_initialization_parameters['simulation']['process_corner']='tt'
	circuit_initialization_parameters['simulation']['conservative']='NO'
	circuit_initialization_parameters['simulation']['w_finger_max']=2e-6

	circuit_initialization_parameters['simulation']['parameters_list']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}

#---------------------------------------------------------------------------------------------------------------------------
# Function that sets the simulation conditions to the optimization_input_parameters dictionary
def get_simulation_conditions_CS_LNA(circuit_initialization_parameters,fo):
	
	circuit_initialization_parameters['simulation']={}

	circuit_initialization_parameters['simulation']['standard_parameters']={}

	# Filenames
	circuit_initialization_parameters['simulation']['standard_parameters']['directory']='/home/ee18b028/cadence_project/lna2/'
	circuit_initialization_parameters['simulation']['standard_parameters']['tcsh']='/home/ee18b028/Optimization/Codes/AutomaticCircuitSynthesis/'
	circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']='series' # 'ideal', 'series','mos_resistor'
	circuit_initialization_parameters['simulation']['standard_parameters']['basic_circuit']='basic_parameters'
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_circuit']='iip3_hb'
	
	# IIP3 Points
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']='basic'		# 'basic' or 'advanced' 
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=-65
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_start']=-70
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_stop']=-40
	circuit_initialization_parameters['simulation']['standard_parameters']['pin_points']=6
	circuit_initialization_parameters['simulation']['standard_parameters']['iip3_calc_points']=3

	# Operating frequency points
	circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=fo
	circuit_initialization_parameters['simulation']['standard_parameters']['f_range']=0

	# Other Values
	circuit_initialization_parameters['simulation']['standard_parameters']['std_temp']=27
	circuit_initialization_parameters['simulation']['standard_parameters']['process_corner']='tt'
	circuit_initialization_parameters['simulation']['standard_parameters']['conservative']='NO'
	circuit_initialization_parameters['simulation']['standard_parameters']['w_finger_max']=2e-6

	circuit_initialization_parameters['simulation']['netlist_parameters']={
		'pin':-65,
		'fund_2':fo+1e6,
		'fund_1':fo,
		'cir_temp':27,
		'n_harm':5
	}


"""
===========================================================================================================================
------------------------------------ Frequency Sweep Code -----------------------------------------------------------------
"""
def frequency_sweep(cir,circuit_parameters,sweep_type,f_start,f_end,n_points,file_location):
	
	# Creating a folder for the output
	if not os.path.exists(file_location):
		os.makedirs(file_location)
		
	# Getting the values of frequency array for the sweep
	if sweep_type=='linear':
		freq_array=np.linspace(f_start,f_end,n_points)
	else:
		freq_array=np.logspace(np.log10(f_start),np.log10(f_end),n_points)

	# Initializing the lists for the analysis
	s11_array=[]
	gain_array=[]
	nf_array=[]
	iip3_array=[]
	Zin_R_array=[]
	Zin_I_array=[]
	vgs_array=[]

	# Starting the sweep
	for freq in freq_array:
		cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=freq
		cir.circuit_initialization_parameters['simulation']['netlist_parameters']['fund_1']=freq
		cir.circuit_initialization_parameters['simulation']['netlist_parameters']['fund_2']=freq+1e6
		cir.update_circuit(circuit_parameters)
		
		
		s11_array.append(cir.extracted_parameters['s11_db'])
		gain_array.append(cir.extracted_parameters['gain_db'])
		nf_array.append(cir.extracted_parameters['nf_db'])
		iip3_array.append(cir.extracted_parameters['iip3_dbm'])
		Zin_R_array.append(cir.extracted_parameters['Zin_R'])
		Zin_I_array.append(cir.extracted_parameters['Zin_I'])
		vgs_array.append(cir.extracted_parameters['vgs_ac'])

	combined_outputs={
		's11':s11_array,
		'Gain':gain_array,
		'NF':nf_array,
		'IIP3':iip3_array,
		'Zin Real':Zin_R_array,
		'Zin Img':Zin_I_array,
		'Vgs':vgs_array
	}

	# Plotting the values for linear sweep
	if sweep_type=='linear':
		for name in combined_outputs:
			figure()
			plot(freq_array,combined_outputs[name])
			grid()
			xlabel('Frequency')
			ylabel(name)
			savefig(file_location+name+'.pdf')
			close()
		
		figure()
		for name in combined_outputs:
			plot(freq_array,combined_outputs[name],label=name)
		grid()
		xlabel('Frequency')
		ylabel('Outputs')
		legend()
		savefig(file_location+'all.pdf')
		close()
	
	# Plotting the values for log sweep
	else:
		for name in combined_outputs:
			figure()
			semilogx(freq_array,combined_outputs[name])
			grid()
			xlabel('Frequency')
			ylabel(name)
			savefig(file_location+name+'.pdf')
			close()
		
		figure()
		for name in combined_outputs:
			semilogx(freq_array,combined_outputs[name],label=name)
		grid()
		xlabel('Frequency')
		ylabel('Outputs')
		legend()
		savefig(file_location+'all.pdf')
		close()

"""
===========================================================================================================================
------------------------------------ Main Program Code --------------------------------------------------------------------
"""

# Creating a dictionary with the optimization parameters
circuit_initialization_parameters={}

# ---------- MOSFET Parameters ----------
#get_mos_parameters(circuit_initialization_parameters,'TSMC180')
get_mos_parameters(circuit_initialization_parameters,'TSMC65')
#get_mos_parameters(circuit_initialization_parameters,'IBM130')

# ---------- Simulation Conditions ----------
fo=1e9
get_simulation_conditions_CS_LNA(circuit_initialization_parameters,fo)

circuit_parameters={
	'Cd':243e-15,
	'Ld':14.2e-9,
	'W':2330e-6,
	'Cg':46.6e-12,
	'Io':303e-6,
	'R1':491,
	'R2':5550,
	'Ls':2.36e-9,
	'Lg':13.1e-9,
	'Rb':5000,
	'Cs':3180e-12,
	
}

f_directory='/home/ee18b028/Optimization/Simulation_Results/CS_LNA/Frequency_Sweep/'
f_name='Test_High_CS_Series_1/'
file_location=f_directory+f_name

sweep_type='linear'
f_start=600e6
f_end=2000e6
n_points=101

cir=sp.Circuit(circuit_initialization_parameters)
frequency_sweep(cir,circuit_parameters,sweep_type,f_start,f_end,n_points,file_location)

#===========================================================================================================================
