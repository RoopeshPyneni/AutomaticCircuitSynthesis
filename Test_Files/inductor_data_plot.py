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
---------------------------------------------- PLOTS -----------------------------------------------------------------------
"""

def plot_single_data(data,circuit_parameters,param_name):
	
	Q_array=[]
	L_array=[]
	param_array=[]

	for i in range(data.shape[0]):
		flag=0
		for param in circuit_parameters:
			if param==param_name:
				continue
			if data.at[i,param]!=circuit_parameters[param]:
				flag=1
				break
		if flag==1:
			continue
		Q_array.append(data.at[i,'Q'])
		L_array.append(data.at[i,'L'])
		param_array.append(data.at[param_name,'Q'])

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
data=pd.read_csv(file_name)
n_rows=data.shape[0]

circuit_parameters={
	'Width':10e-6,
	'Radius':30e-6,
	'N_Turns':1.25,
	'gdis':20e-6,
	'spc':3e-6
}
# f.write('Width,Radius,N_Turns,gdis,spc,Q,L\n')

plot_single_data(data,circuit_parameters,'Width')


