#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform the iip3 analysis by varying pin and plotting the fundamental and im3 component of vout
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf
from matplotlib import pylab
from pylab import *
#===========================================================================================================================


#===========================================================================================================================
#------------------------------------Defining the functions -----------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit_parameters to a txt file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def write_circuit_parameters(circuit_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']	# Getting the filename
	
	newpath =filename+'/Sensitivity_Analysis/Results'	# Creating the folder if it is not present
	if not os.path.exists(newpath):
		os.makedirs(newpath)

	filename=filename+'/Sensitivity_Analysis/Results/circuit_parameters.txt'
	
	f=open(filename,'w')
	for param_name in circuit_parameters:
		f.write(str(param_name)+'\t'+str(circuit_parameters[param_name])+'\n')	# Writing the values in the file
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for extracted parameters to a csv file
# Input: extracted_parameters, optimization_input_parameters
# Output: NONE
def write_extracted_parameters_initial(optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/iip3.csv'	# Getting the filename

	f=open(filename,'w')
	f.write('Frequency,Pin,Vout_Fund,Vout_IM3,IIP3,Fund_Slope,IM3_Slope\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_extracted_parameters(freq,extracted_parameters,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/iip3.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(freq))
	f.write(','+str(extracted_parameters['iip3_pin']))
	f.write(','+str(extracted_parameters['iip3_fund']))
	f.write(','+str(extracted_parameters['iip3_im3']))
	f.write(',,,')
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of extracted_parameters from each temperature iteration to a csv file
# Input: circuit_parameters, optimization_input_parameters
# Output: NONE
def update_iip3(freq,iip3,fund_slope,im3_slope,optimization_input_parameters):
	
	filename=optimization_input_parameters['filename']['output']+'/Sensitivity_Analysis/Results/iip3.csv'	# Getting the filename
	
	f=open(filename,'a')
	f.write(str(freq))
	f.write(',')
	f.write(',')
	f.write(',')
	f.write(','+str(iip3)+','+str(fund_slope)+','+str(im3_slope))
	f.write('\n')
	f.close()
	


"""
===========================================================================================================================
--------------------------------------------Output Functions---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function that will perform the iip3 analysis
def iip3_analysis(cir,optimization_input_parameters,timing_results):
	
	if optimization_input_parameters['iip3_analysis']['run']=='NO':
		return

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('IIP3 Analysis Start\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()

	# Storing the starting time
	timing_results['iip3_analysis']={}
	timing_results['iip3_analysis']['start']=datetime.datetime.now()

	print('************************************************************************************************************')
	print('****************************************** IIP3 Analysis ***************************************************')
	
	# Storing the initial values
	initial_freq=cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']
	initial_pin=cir.circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']
	cir.circuit_initialization_parameters['simulation']['standard_parameters']['iip3_type']='basic'
	cir.update_simulation_parameters(optimization_input_parameters['iip3_analysis']['simulation'])
	
	# Pin values and frequency array
	pin_start=optimization_input_parameters['iip3_analysis']['pin_start']
	pin_stop=optimization_input_parameters['iip3_analysis']['pin_stop']
	n_pin=optimization_input_parameters['iip3_analysis']['n_pin']
	n_points=optimization_input_parameters['iip3_analysis']['n_points']
	freq_array=optimization_input_parameters['iip3_analysis']['freq_array']

	pin_array=np.linspace(pin_start,pin_stop,n_pin)
	
	# Writing the values to output files
	write_circuit_parameters(cir.circuit_parameters,optimization_input_parameters)
	write_extracted_parameters_initial(optimization_input_parameters)
	
	# Performing the analysis
	for freq in freq_array:
		cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=freq
		im3_array=[]
		fund_array=[]
		for pin in pin_array:
			cir.circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=pin
			cir.run_circuit()
			fund_array.append(cir.extracted_parameters['iip3_fund'])
			im3_array.append(cir.extracted_parameters['iip3_im3'])
			update_extracted_parameters(freq,cir.extracted_parameters,optimization_input_parameters)		# Writing the values to the output file
	
		# Updating the calculated iip3 slope
		fund_array=np.array(fund_array)
		im3_array=np.array(im3_array)
		iip3,_,im3_slope,_,fund_slope=cir.calculate_iip3(n_pin,n_points,cf.dbm_to_normal(fund_array),cf.dbm_to_normal(im3_array),pin_array)
		update_iip3(freq,iip3,fund_slope,im3_slope,optimization_input_parameters)

		# Plotting the graphs
		file_directory=optimization_input_parameters['filename']['output']
		plot_iip3_analysis(freq,pin_array,im3_array,fund_array,file_directory)
		
	# Restoring the value of initial extracted and circuit parameters
	cir.circuit_initialization_parameters['simulation']['standard_parameters']['pin_fixed']=initial_pin
	cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_operating']=initial_freq
	cir.run_circuit()

	# Storing the starting time
	timing_results['iip3_analysis']['stop']=datetime.datetime.now()

	# Opening the Run_Status File
	f=open(optimization_input_parameters['filename']['run_status'],'a')
	f.write('Sensitivity Analysis End\n Time : '+str(datetime.datetime.now())+'\n\n')
	f.close()
	
	
#===========================================================================================================================


"""
====================================================================================================================================================================
"""

#-----------------------------------------------------------------------------------------------
# Plotting results ( Parameters vs Io at different temperatures )
# Inputs  : extracted_parameters_iter
# Outputs : extracted_matrix, temp_array, current_array, param_array
def plot_iip3_analysis(freq,pin_array,im3_array,fund_array,file_directory):

	# Creating a directory
	file_directory=file_directory+'/IIP3_Analysis/Plots/'
	if not os.path.exists(file_directory):
		os.makedirs(file_directory)
	
	figure()
	plot(pin_array,fund_array,color='green',label='Fundamental')
	plot(pin_array,im3_array,color='red',label='IM3')
	xlabel('Pin')
	ylabel('Pout')
	legend()
	grid()
	savefig(file_directory+str(freq)+'.pdf')
	close()


#===========================================================================================================================
