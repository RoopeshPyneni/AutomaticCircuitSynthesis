#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file is used to plot the values of Q and L from the csv file where data is stored
"""

#===========================================================================================================================
import numpy as np
import fileinput
import os
import pandas as pd
from matplotlib import pylab
from pylab import *


"""
============================================================================================================================
---------------------------------------------- DATA EXTRACTION -------------------------------------------------------------
"""

def extract_data(file_name):

	# Reading the data from the csv file
	f=open(file_name,'r')
	lines=f.readlines()
	f.close()

	# Initializing a dictionary
	output_dict={}
	header_dict={}

	# Reading line 0 ( HEADER )
	line=lines[0]
	lines=lines[1:]
	line=line.split('\n')[0]
	line=line.split(',')
	i=0
	for name in line:
		output_dict[name]=[]
		header_dict[i]=name
		i+=1

	# Reading other lines
	n_rows=0
	while lines!=[]:
		n_rows+=1
		line=lines[0]
		lines=lines[1:]
		line=line.split(',')
		i=0
		for val in line:
			output_dict[header_dict[i]].append(float(val))
			i+=1
	
	return output_dict,n_rows


"""
============================================================================================================================
---------------------------------------------- PLOTS -----------------------------------------------------------------------
"""

def plot_single_data(data,n_rows,circuit_parameters,param_name):
	
	Q_array=[]
	L_array=[]
	param_array=[]

	for i in range(n_rows):
		flag=0
		for param in circuit_parameters:
			if param==param_name:
				continue
			if data[param][i]!=circuit_parameters[param]:
				flag=1
				break
		if flag==1:
			continue
		Q_array.append(data['Q'][i])
		L_array.append(data['L'][i])
		param_array.append(data[param_name][i])

	figure()
	plot(param_array,L_array)
	xlabel('Parameter Value')
	ylabel('L')
	grid()
	show()
	close()

	figure()
	plot(param_array,Q_array)
	xlabel('Parameter Value')
	ylabel('Q')
	grid()
	show()
	close()



"""
============================================================================================================================
---------------------------------------------- MAIN PROGRAM ----------------------------------------------------------------
"""

# Filenames for csv file
file_directory_output='/home/ee18b028/Optimization/Simulation_Results/Inductor/Sweep/'
file_name=file_directory_output+'inductor_sweep.csv'

# Reading the csv data
data,n_rows=extract_data(file_name)

circuit_parameters={
	'Width':10e-6,
	'Radius':30e-6,
	'N_Turns':1.25,
	'gdis':20e-6,
	'spc':3e-6
}
# f.write('Width,Radius,N_Turns,gdis,spc,Q,L\n')

plot_single_data(data,circuit_parameters,'Width')


