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
------------------------------------------ OTHER FUNCTION ------------------------------------------------------------------
"""

#--------------------------------------------------------------------------------
# This function will sort two lists
def sort_2(list1,list2,**kwargs):
    list_zip=list(zip(list1,list2))
    if 'order' in kwargs:
        if kwargs['order']=='D':
            list_zip=sorted(list_zip,reverse=True)
        else:
            list_zip=sorted(list_zip,reverse=False)
    else:
        list_zip=sorted(list_zip,reverse=False)
    list1=[i for (i,s) in list_zip]
    list2=[s for (i,s) in list_zip]
    return list1,list2
    
    
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
		column_dict[n_col]=name
		data_dict[name]=[]
		n_col+=1

	# Getting the data
	n_rows=0
	while lines!=[]:
		line=lines[0]
		lines=lines[1:]
		line=line.split('\n')[0]
		line=line.split(',')
		i=0
		for name in line:
			data_dict[column_dict[i]].append(float(name))
			i+=1
		n_rows+=1
		
		if n_rows%10000==0:
			print(n_rows)
	
	#print(data_dict)
	print(n_rows)

	return data_dict

#---------------------------------------------------------------------------------------------------------------------------
# Sorting the data
def sort_data(data_dict):
	
	L_array=data_dict['L'].copy()
	new_L=np.arange(len(L_array))
	
	L_array,new_L=sort_2(L_array,new_L)
	
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
file_directory='/home/ee18b028/Optimization/Simulation_Results/Inductor/Sweep2/'

print('Starting Data Extract')
print(datetime.datetime.now())
data_dict=extract_data(file_directory+'inductor_sweep.csv')

print('Starting Data Sort')
print(datetime.datetime.now())
new_L=sort_data(data_dict)

print('Starting Data Write')
print(datetime.datetime.now())
store_data(file_directory+'inductor_sweep_1.csv',data_dict,new_L)






