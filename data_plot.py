#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Data Plot File:
"""
#===========================================================================================================================
import numpy as np
import os
import common_functions as cf
from matplotlib import pylab
from pylab import *
#===========================================================================================================================


	

#===========================================================================================================================
#----------------------------------- Extracting and Plotting Double Array --------------------------------------------------
	
#-----------------------------------------------------------------------------------------------
# Function to extract the loss data from csv file
def extract_double_array(filename_root,filename_name):
	
	# Extracting data from the files
	filename=filename_root+'/results/'+filename_name+'.csv'
	f=open(filename)
	lines=f.readlines()
	f.close()
	
	# Assigning the array to store the loss variables
	param_name=[]
	variable_name=[]
	
	# Extracting the first line of the csv file
	line=lines[0].split(',')
	lines=lines[1:]
	
	# Creating variables for the counting no of lines and variables
	n_iter=len(lines)
	n_param=0
	n_variable=0
	
	# Counting n_param and n_variable
	for char in line:
		if char == '\n':
			break
		n_param+=1
		
	n_param-=1 # CHANGE
	
	for i in range(3,n_param):
		if line[i]==line[2]:
			n_variable=i-3
			break

	n_param=int(n_param/n_variable)
	
	
	# Storing n_param and n_variable
	for i in range(n_variable):
		variable_name.append(line[i+2])
	for i in range(n_param):
		param_name.append(line[(n_variable+1)*i+1])	
	
	# Creating array to store the values 
	loss_slope_array=np.zeros((n_iter,n_param,n_variable),dtype=float)
	
	# Extracting the values of the variables
	for i in range(n_iter):
	
		# Extracting the next line
		line=lines[0].split(',')
		lines=lines[1:]
		
		# Storing the variables
		for j in range(n_param):
			for k in range(n_variable):
				loss_slope_array[i,j,k]=float(line[(n_variable+1)*j+k+2])
		
	return loss_slope_array,param_name,variable_name
	
#-----------------------------------------------------------------------------------------------
# Plotting loss slope vs iterations
def plot_double_array(filename_root,filename_name):
	
	loss_array,param_name,variable_name=extract_double_array(filename_root,filename_name)
	n_iter=loss_array.shape[0]
	n_param=loss_array.shape[1]
	n_variable=loss_array.shape[2]
	
	filename=filename_root+'/plots/'+filename_name+'/'
	if not os.path.exists(filename):
    		os.makedirs(filename)
	
	# Creating the new arrays
	arrX=np.zeros((n_iter,1),dtype=float)
	
	# Calculating values of new array
	i=0	
	while i<n_iter:
		arrX[i,0]=i+1	
		i+=1
		
	print('Starting Plots '+filename_name)
	
	color_dict={0:'r',1:'g',2:'b',3:'c',4:'m',5:'k',6:'y'}
	n_colour=7
	
	for i in range(n_param):
		figure()
		for j in range(n_variable):
			plot(arrX,loss_array[:,i,j],color_dict[int(j%n_colour)],label=variable_name[j])
		xlabel('Iterations')
		ylabel(param_name[i])
		legend()
		grid()
		savefig(filename+param_name[i]+'.pdf')
		close()
	
	print('Plotting Over '+filename_name)

	

#===========================================================================================================================
#------------------------------------Extracting and Plotting Single Array --------------------------------------------------

#-----------------------------------------------------------------------------------------------
# Function to extract the data from csv file
def extract_single_array(filename_root,filename_name):
	
	# Extracting data from the files
	filename=filename_root+'/results/'+filename_name+'.csv'
	f=open(filename)
	lines=f.readlines()
	f.close()
	
	# Assigning the array to store the loss variables
	param_name=[]
	
	# Extracting the first line of the csv file
	line=lines[0].split(',')
	lines=lines[1:]
	
	# Creating variables for the counting no of lines and variables
	n_iter=len(lines)
	n_param=0
	
	# Storing the variable names
	for char in line:
		if char == '\n':
			break
		if char != 'Iteration No':
			param_name.append(char)
		n_param+=1

	n_param-=1 # CHANGE
		
	# Creating array to store the values 
	loss_array=np.zeros((n_iter,n_param),dtype=float)
	
	# Extracting the values of the variables
	for i in range(n_iter):
	
		# Extracting the next line
		line=lines[0].split(',')
		lines=lines[1:]
		
		# Storing the variables
		for j in range(n_param):
			loss_array[i,j]=float(line[j+1])
		
	return loss_array,param_name
	
	
#-----------------------------------------------------------------------------------------------
# Plotting average parameters vs iterations
def plot_single_array(filename_root,filename_name):
	
	loss_array,param_name=extract_single_array(filename_root,filename_name)
	n_iter=loss_array.shape[0]
	n_param=loss_array.shape[1]
	
	filename=filename_root+'/plots/'+filename_name+'/'
	if not os.path.exists(filename):
    		os.makedirs(filename)
	
	# Creating the new arrays
	arrX=np.zeros((n_iter,1),dtype=float)
	
	# Calculating values of new array
	i=0	
	while i<n_iter:
		arrX[i,0]=i+1	
		i+=1
		
	print('Starting Plots '+filename_name)
	
	color_dict={0:'r',1:'g',2:'b',3:'c',4:'m',5:'k',6:'y'}
	n_colour=7
	
	# Figure 1
	figure()
	for i in range(n_param):
		plot(arrX,loss_array[:,i],color_dict[int(i%n_colour)],label=param_name[i])
		annotate(cf.num_trunc(loss_array[n_iter-1,i],3),(arrX[n_iter-1,0],loss_array[n_iter-1,i]))
	xlabel('Iterations')
	ylabel('Parameter')
	legend()
	grid()
	savefig(filename+'all.pdf')
	close()
	
	
	# Figure 2	
	for i in range(n_param):
		figure()
		plot(arrX,loss_array[:,i],'g',label=param_name[i])
		annotate(cf.num_trunc(loss_array[n_iter-1,i],3),(arrX[n_iter-1,0],loss_array[n_iter-1,i]))
		xlabel('Iterations')
		ylabel(param_name[i])
		legend()
		grid()
		savefig(filename+param_name[i]+'.pdf')
		close()
	
	print('Plotting Over '+filename_name)

#-----------------------------------------------------------------------------------------------
# Function to extract the data from temperature analysis csv file
def extract_temp_analysis(filename):
	
	# Extracting data from the files
	f=open(filename)
	lines=f.readlines()
	f.close()
	
	# Assigning the array to store the loss variables
	param_name=[]
	
	# Extracting the first line of the csv file
	line=lines[0].split(',')
	lines=lines[1:]
	
	# Creating variables for the counting no of lines and variables
	n_iter=len(lines)
	n_param=0
	
	# Storing the variable names
	for char in line:
		if char == '\n':
			break
		param_name.append(char)
		n_param+=1
	
	# Creating array to store the values 
	extracted_array=np.zeros((n_iter,n_param),dtype=float)
	
	# Extracting the values of the variables
	for i in range(n_iter):
	
		# Extracting the next line
		line=lines[0].split(',')
		lines=lines[1:]
		
		# Storing the variables
		for j in range(n_param):
			extracted_array[i,j]=float(line[j])
		
	return extracted_array,param_name

#-----------------------------------------------------------------------------------------------
# Plotting results from temperature analysis
def plot_temp_analysis(filename):
	extract_filename=filename+'/Temperature_Analysis/Results/extracted_parameters.csv'
	extracted_array,param_name=extract_temp_analysis(extract_filename)
	pathname=filename+'/Temperature_Analysis/Plots/'
	filename=filename+'/Temperature_Analysis/Plots/'
	
	if not os.path.exists(pathname):
    		os.makedirs(pathname)
    		
	n_param=len(param_name)

	arrX=extracted_array[:,0]

	output_conditions_dict={'gain_db':10.0,'s11_db':-15.0,'nf_db':4.0,'iip3_dbm':-5.0}

	# Figure 1
	for i in range(n_param-1):
		figure()
		plot(arrX,extracted_array[:,i+1],'g',label=param_name[i+1])
		if param_name[i+1] in output_conditions_dict:
			arrY=np.ones(len(arrX),dtype=float)
			arrY=arrY*output_conditions_dict[param_name[i+1]]
			plot(arrX,arrY,'r',label='Output Condition')
		xlabel('Temperature')
		ylabel(param_name[i+1])
		legend()
		grid()
		savefig(filename+param_name[i+1]+'.pdf')
		close()



	

#===========================================================================================================================
#------------------------------------Main Program Code----------------------------------------------------------------------
def plot_complete(optimization_input_parameters):

	filename_temp=optimization_input_parameters['filename']['output']
	filename=filename_temp+'/Optimization'
	
	plot_single_array(filename,'loss')
	plot_single_array(filename,'alpha_parameters')
	plot_single_array(filename,'output_parameters')
	plot_single_array(filename,'circuit_parameters')
	plot_single_array(filename,'average_parameters')
	
	plot_double_array(filename,'loss_slope')
	plot_double_array(filename,'sensitivity')

	plot_temp_analysis(filename_temp)
	

#===========================================================================================================================

