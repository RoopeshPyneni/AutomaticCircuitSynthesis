#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file is used to store the value of Q and L of TSMC Inductor by sweeping various inductor parameters
"""

#===========================================================================================================================
import numpy as np
import fileinput
import os
from matplotlib import pylab
from pylab import *


"""
============================================================================================================================
------------------------------------------ EXTRACTION FUNCTION -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Extracting the data from the file
def extract_data(file_name):
	f=open(file_name)
	lines=f.readlines()
	f.close()

	# Variables
	data_dict={}
	column_dict={}
	n_col=0

	# Getting the field names
	line=lines[0]
	lines=lines[1:]
	line=line.split('\n')[0]
	line=line.split(',')
	for name in line:
		column_dict[n_col]=column_dict
		data_dict[name]=[]
		n_col+=1

	# Getting the data
	while lines!=[]:
		line=lines[0]
		lines=lines[1:]
		line=line.split('\n')[0]
		line=line.split(',')
		i=0
		for name in line:
			data_dict[column_dict[i]].append(float(name))
			i+=1

	return data_dict

#---------------------------------------------------------------------------------------------------------------------------
# Sorting the data
def sort_data(data_dict):
	L_array=data_dict['L']
	n_rows=len(L_array)

	new_L=[]
	for i in range(n_rows):
		i_min=-1
		L_min=1e15
		for j in range(n_rows):
			if j in new_L:
				continue
			if L_array[j]<L_min:
				i_min=j
				L_min=L_array[j]
		
		if i_min>-1:
			new_L.append(i_min)
	
	return new_L

#---------------------------------------------------------------------------------------------------------------------------
# Storing the data
def store_data(filename,data_dict,new_L):

	f=open(filename,'w')

	# Writing the column names
	i=0
	for key in data_dict:
		if i>0:
			f.write(',')
		elif i==0:
			i=1
		f.write(str(key))
	f.write('\n')

	# Writing the data
	for i in range(len(new_L)):
		flag=0
		for key in data_dict:
			if flag>0:
				f.write(',')
			elif flag==0:
				flag=1
			f.write(str(data_dict[key][new_L[i]]))
		f.write('\n')
	
	f.close()



"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_directory='/home/ee18b028/cadence_project/Test_Circuits/Inductor/Inductor_Test_1/'
data_dict=extract_data(file_directory+'inductor_sweep.csv')
new_L=sort_data(data_dict)
store_data(file_directory+'inductor_sweep_1.csv',data_dict,new_L)






