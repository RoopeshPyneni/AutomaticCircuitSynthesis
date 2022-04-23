#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file is used to store the value of Q and L of TSMC Inductor by sweeping various inductor parameters
"""

#===========================================================================================================================
from operator import length_hint
import numpy as np
import pandas as pd
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
# Sorting the data
def sort_data(data):
	
	width=data['Width'].tolist()
	length=data['Length'].tolist()
	capacitance=data['Capacitance'].tolist()
	
	new_cap=np.arange(len(capacitance))
	capacitance,new_cap=sort_2(capacitance,new_cap)

	width=[width[i] for i in new_cap]
	length=[length[i] for i in new_cap]

	data['Width']=width
	data['Length']=length
	data['Capacitance']=capacitance
	
	return data


"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

# Filenames for the netlist file
file_name='/home/ee18b028/Optimization/Simulation_Results/Capacitor/Sweep2/capacitor_sweep.csv'


print('Starting Data Extract')
print(datetime.datetime.now())
data=pd.read_csv(file_name)

print('Starting Data Sort')
print(datetime.datetime.now())
data=sort_data(data)

print('Starting Data Write')
print(datetime.datetime.now())
data.to_csv(file_name)






