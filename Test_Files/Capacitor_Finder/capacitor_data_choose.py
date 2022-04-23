#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file is used to store the value of Q and L of TSMC Inductor by sweeping various inductor parameters
"""

#===========================================================================================================================
import numpy as np
from matplotlib import pylab
from pylab import *
import pandas as pd


"""
============================================================================================================================
------------------------------------------ FUNCTIONS -----------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Finding the location of L
def find_target_L(data,L,start,end):

	n1=int((start+end+1)/2)
	
	if data.at[n1,'L']>=L and data.at[n1-1,'L']<=L:
		return n1
	elif data.at[n1,'L']>=L and data.at[n1-1,'L']>=L:
		return find_target_L(data,L,start,n1-1)
	else:
		return find_target_L(data,L,n1,end)


#---------------------------------------------------------------------------------------------------------------------------
# Finding the best point
def find_best_points(Q,L):

	# Reading the data
	file_directory='/home/ee18b028/Optimization/Simulation_Results/Inductor/Sweep2/'
	#file_directory='C:/Users/roope/Studies/IIT/Prof Projects/Circuit_Synthesis/Extra_Codes/'
	data=pd.read_csv(file_directory+'inductor_sweep_1.csv')
	n_rows=data.shape[0]

	# Finding the location of 0.98L and 1.02L
	print(n_rows-1)
	n1=find_target_L(data,L*0.98,0,n_rows-1)
	n2=find_target_L(data,L*1.02,0,n_rows-1)

	# Finding the loss of all the points from n1 and n2
	loss_iter=[]
	for i in range(n1,1+n2):
		Qi=data.at[i,'Q']
		Li=data.at[i,'L']
		loss=((Qi-Q)/Q)**2
		loss+=((Li-L)/L)**2
		loss_iter.append(loss)

	# Finding the location of the smallest value in loss iter
	min_value=min(loss_iter)
	min_index=loss_iter.index(min_value)
	min_index+=n1

	return min_index,data.at[min_index,'Width'],data.at[min_index,'Radius'],data.at[min_index,'N_Turns'],data.at[min_index,'gdis'],data.at[min_index,'spc']


"""
====================================================================================================================================================================================
------------------------------------------------------------ MAIN PROGRAM ----------------------------------------------------------------------------------------------------------
"""

Q=10
L=1e-9
print(find_best_points(Q,L))



