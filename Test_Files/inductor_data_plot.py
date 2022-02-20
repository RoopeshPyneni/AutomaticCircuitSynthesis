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
from matplotlib import pylab
from pylab import *

"""
============================================================================================================================
---------------------------------------------- OTHER FUNCTIONS -------------------------------------------------------------
"""

def get_unique(arr):
	arr1=[]
	for a in arr:
		if a not in arr1:
			arr1.append(a)
	return arr1

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
	
	# Printing the unique values of each field
	print(len(output_dict['Width']))
	print(n_rows)
	print(get_unique(output_dict['Width']))
	print(get_unique(output_dict['Radius']))
	print(get_unique(output_dict['N_Turns']))
	print(get_unique(output_dict['gdis']))
	print(get_unique(output_dict['spc']))
		
	return output_dict,n_rows


"""
============================================================================================================================
---------------------------------------------- PLOTS -----------------------------------------------------------------------
"""

def plot_single_data(data,n_rows,circuit_parameters,param_name,file_directory):
	
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
	
	print(Q_array)
	print(L_array)
	print(param_array)

	figure()
	plot(param_array,L_array)
	xlabel('Parameter Value')
	ylabel('L')
	grid()
	#show()
	savefig(file_directory+'L.pdf')
	close()

	figure()
	plot(param_array,Q_array)
	xlabel('Parameter Value')
	ylabel('Q')
	grid()
	#show()
	savefig(file_directory+'Q.pdf')
	close()

def plot_scatter_Q_L(data,n_rows,file_directory):
	
	Q_array=[]
	L_array=[]

	for i in range(n_rows):
		Q_array.append(data['Q'][i])
		L_array.append(data['L'][i])
		
	figure()
	scatter(Q_array,np.log10(L_array),s=2)
	#plot(param_array,L_array)
	xlabel('Q Value')
	ylabel('L Value')
	grid()
	#show()
	savefig(file_directory+'Q_L.pdf')
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
	'Width':9e-6,
	'Radius':40e-6,
	'N_Turns':1.5,
	'gdis':40e-6,
	'spc':3e-6
}
# f.write('Width,Radius,N_Turns,gdis,spc,Q,L\n')

print('Minimum L :',min(data['L']))
print('Maximum L :',max(data['L']))
print('Minimum Q :',min(data['Q']))
print('Maximum Q :',max(data['Q']))

plot_single_data(data,n_rows,circuit_parameters,'spc',file_directory_output)
#plot_scatter_Q_L(data,n_rows,file_directory_output)


