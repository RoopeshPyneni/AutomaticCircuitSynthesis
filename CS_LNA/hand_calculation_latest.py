#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform hand calculations for CS LNA
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf # type: ignore
import spectre_common as sc # type: ignore


"""
===========================================================================================================================
------------------------------------Defining the functions for simple calculations-----------------------------------------
"""


#-----------------------------------------------------------------------------------------------
# Converting from dB to normal
def db_to_normal(val_db):
    return 10**(val_db/10)

#-----------------------------------------------------------------------------------------------
# Calculating Cd
def calculate_Cd(Cload,Ld,fo):
	wo=2*np.pi*fo
	Cd=1/(wo*wo*(Ld))
	Cd-=Cload
	return Cd

#-----------------------------------------------------------------------------------------------
# Calculating Qin
def calculate_Qin(s11,fo,f_delta):
	Ks=10**(-0.2*s11)-1
	
	f1=1+(f_delta/fo)
	f2=1-(f_delta/fo)

	kf1=f1-(1/f1)
	kf2=f2-(1/f2)

	kf=max(kf1,kf2)

	Qin=(1/kf)*np.sqrt(1/Ks)

	return Qin

#-----------------------------------------------------------------------------------------------
# Calculating Cgs
def calculate_cgs(Rs,fo,Qin):
    wo=2*np.pi*fo
    cgs=1/(2*Rs*Qin*wo)
    return cgs

#-----------------------------------------------------------------------------------------------
# Calculating gm
def calculate_gm(Ld,fo,cgs,nf):
	wo=2*np.pi*fo
	Rd=wo*Ld*15
	A=200*((wo*cgs)**2)/Rd
	B=100*((wo*cgs)**2)
	C=1-10**(0.1*nf)
	return 8*A/(np.sqrt(B*B-4*A*C)-B)


#-----------------------------------------------------------------------------------------------
# Calculating W
def calculate_W(cgs,Lmin,Cox):
    return 3*cgs/(2*Lmin*Cox)

#-----------------------------------------------------------------------------------------------
# Calculating Io
def calculate_Io(gm,un,cox,W,Lmin):
    return (gm*gm)/(2*un*cox*W/Lmin)

#-----------------------------------------------------------------------------------------------
# Calculating Ls
def calculate_Ls(rsource,cgs,gm,fo):
    return rsource*cgs/gm
    
#-----------------------------------------------------------------------------------------------
# Calculating Lg
def calculate_Lg(Ls,cgs,fo):
    w=2*np.pi*fo
    Lg=1/(w*w*cgs)-Ls
    return Lg

#-----------------------------------------------------------------------------------------------
# Calculating Zim_max
def calculate_Zim_max(s11):
	s11_mag_sq=10**(0.1*s11)
	H_by_Zim_max_sq=(1/s11_mag_sq)-1
	Z_im_max=100*np.sqrt(1/H_by_Zim_max_sq)
	return Z_im_max

#-----------------------------------------------------------------------------------------------
# Calculating R1 and R2
def calculate_Rsum_Rk(vdd,fo,Cg,Imax):

	# Getting the value of k
	Rk=0.9 # This is set arbitrarily

	# Getting the value of Rsum from Reff
	Reff=10/(2*np.pi*fo*Cg)
	Rsum=Reff/(Rk*(1-Rk))
	
	# Multiplying R1 and R2 to get a higher value
	Rmin=vdd/Imax
	if Rsum<Rmin:
		Rsum=Rmin

	return Rsum,Rk

#-----------------------------------------------------------------------------------------------
# Updating R1 and R2
def updating_R1_R2(vdsat1,vdsat2,vth,vdd,Imax):
	vg=vdsat1+vdsat2+vth
	R2=1e4*vg/vdd
	R1=1e4-R2

	# Multiplying R1 and R2 to get a higher value
	Rmin=vdd/Imax
	if R1+R2<Rmin:
		k=Rmin/(R1+R2)
		R1*=k
		R2*=k

	return R1,R2

#-----------------------------------------------------------------------------------------------
# Updating Cd
def updating_Cd(Cgd,Cload,Ld,fo):
	wo=2*np.pi*fo
	Cd=1/(wo*wo*(Ld))
	Cd-=Cload
	Cd-=Cgd
	return Cd

#-----------------------------------------------------------------------------------------------
# Updating gm
def update_gm(extracted_nf,target_nf,gm):

	if extracted_nf<(target_nf-0.3):
		return gm/1.2
	elif extracted_nf>target_nf:
		return gm*1.2
	else:
		return gm

#---------------------------------------------------------------------------------------------------------------------------
# Getting the best point
def get_best_point(circuit_parameters_iter,extracted_parameters_iter,output_conditions):
	
	# Getting the best results for the final pre-optimized point
	j_array=[key for key in circuit_parameters_iter]
	loss_array=[]
	loss_zero=0
	Io_array=[]
	
	# Getting the value of Io and loss
	for j in j_array:
		Io_array.append(extracted_parameters_iter[j]['Io'])
		loss_gain=sc.ramp_func(output_conditions['gain_db']-extracted_parameters_iter[j]['gain_db'])
		loss_s11=sc.ramp_func(extracted_parameters_iter[j]['s11_db']-output_conditions['s11_db'])
		loss_nf=sc.ramp_func(extracted_parameters_iter[j]['nf_db']-output_conditions['nf_db'])
		loss_array.append(loss_gain+loss_nf+loss_s11)
		if loss_array[-1]==0:
			loss_zero=1
	
	# Checking for iter number if loss=0
	if loss_zero==1:
		Io_array=[Io_array[j] for j in j_array if loss_array[j]==0]
		j_array=[j for j in j_array if loss_array[j]==0]
		return j_array[Io_array.index(min(Io_array))]
	else:
		return j_array[Io_array.index(min(Io_array))]
	


"""
===========================================================================================================================
-------------------------------------------- Storing Results --------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Writing the header row for circuit parameters and extracted parameters to a csv file
def write_parameters_initial(cir,optimization_input_parameters):
	
	# Creating a file path
	filepath=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/'
	if not os.path.exists(filepath):
		os.makedirs(filepath)


	# Storing results for circuit parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/initial_circuit_parameters.csv'

	f=open(filename,'w')
	f.write('Iter_Number,Iter_Type')
	for param_name in cir.get_initial_circuit_parameters():
		f.write(','+param_name)
	f.write('\n')
	f.close()


	# Storing results for circuit parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/circuit_parameters.csv'

	f=open(filename,'w')
	f.write('Iter_Number,Iter_Type')
	for param_name in cir.get_circuit_parameters():
		f.write(','+param_name)
	f.write('\n')
	f.close()


	# Storing results for extracted parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/extracted_parameters.csv'

	f=open(filename,'w')
	f.write('Iter_Number,Iter_Type')
	for param_name in cir.get_extracted_parameters():
		f.write(','+param_name)
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit parameters and extracted parameters from each iteration to a csv file
def update_parameters(cir,optimization_input_parameters,iter_no,iter_type):
	
	# Storing results for circuit parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/initial_circuit_parameters.csv'
	
	f=open(filename,'a')
	f.write(str(iter_no)+','+str(iter_type))
	for param_name in cir.get_initial_circuit_parameters():
		f.write(','+str(cir.get_initial_circuit_parameters()[param_name]))
	f.write('\n')
	f.close()

	# Storing results for circuit parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/circuit_parameters.csv'
	
	f=open(filename,'a')
	f.write(str(iter_no)+','+str(iter_type))
	for param_name in cir.get_circuit_parameters():
		f.write(','+str(cir.get_circuit_parameters()[param_name]))
	f.write('\n')
	f.close()

	# Storing results for extracted parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/extracted_parameters.csv'
	
	f=open(filename,'a')
	f.write(str(iter_no)+','+str(iter_type))
	for param_name in cir.get_extracted_parameters():
		f.write(','+str(cir.get_extracted_parameters()[param_name]))
	f.write('\n')
	f.close()


"""
===========================================================================================================================
-------------------------------------------- Main Functions ---------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the Initial Circuit Parameters
def calculate_initial_parameters(cir,optimization_input_parameters):

	output_conditions=optimization_input_parameters['output_conditions']
    
	# Getting the output conditions
	Cload=output_conditions['Cload']
	fo=output_conditions['wo']/(2*np.pi)
	f_range=cir.circuit_initialization_parameters['simulation']['standard_parameters']['f_range']
	Rs=output_conditions['Rs']
	s11=output_conditions['s11_db']
	nf=output_conditions['nf_db']
	vdd=cir.mos_parameters['vdd']
	Lmin=cir.mos_parameters['Lmin']
	Cox=cir.mos_parameters['cox']
	un=cir.mos_parameters['un']
	vth=cir.mos_parameters['vt']

	# Calculating the circuit parameters
	initial_circuit_parameters={}
	initial_circuit_parameters['Ld']=9e-9
	initial_circuit_parameters['Cd']=calculate_Cd(Cload,initial_circuit_parameters['Ld'],fo)
	Qin=calculate_Qin(s11,fo,f_range)
	cgs=calculate_cgs(Rs,fo,Qin)
	initial_circuit_parameters['W']=calculate_W(cgs,Lmin,Cox)
	global gm
	gm=calculate_gm(initial_circuit_parameters['Ld'],fo,cgs,nf)
	initial_circuit_parameters['Cg']=100*cgs
	initial_circuit_parameters['Io']=calculate_Io(gm,un,Cox,initial_circuit_parameters['W'],Lmin)
	initial_circuit_parameters['Rsum'],initial_circuit_parameters['Rk']=calculate_Rsum_Rk(vdd,fo,initial_circuit_parameters['Cg'],optimization_input_parameters['pre_optimization']['I_Rdivider_max'])
	initial_circuit_parameters['Ls']=calculate_Ls(Rs,cgs,gm,fo)
	initial_circuit_parameters['Lg']=calculate_Lg(initial_circuit_parameters['Ls'],cgs,fo)
	initial_circuit_parameters['Rb']=5000
	initial_circuit_parameters['Cs']=100/(2*np.pi*50*fo)
	if cir.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_capacitor' or cir.circuit_initialization_parameters['simulation']['standard_parameters']['circuit_type']=='mos_inductor':
		initial_circuit_parameters['Cs']=10/(2*np.pi*50*fo)

	# Running the circuit
	cir.update_circuit(initial_circuit_parameters)

#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters	after calculating the new value of vt
def update_initial_parameters(cir,optimization_input_parameters):

	i=0
	Lmin=cir.mos_parameters['Lmin']
	Cox=cir.mos_parameters['cox']
	un=cir.mos_parameters['un']
	vdd=cir.mos_parameters['vdd']
  
    	# Getting the output conditions
	Cload=optimization_input_parameters['output_conditions']['Cload']
	fo=optimization_input_parameters['output_conditions']['wo']/(2*np.pi)
	nf=optimization_input_parameters['output_conditions']['nf_db']

	write_parameters_initial(cir,optimization_input_parameters)

	initial_circuit_parameters_iter={}
	circuit_parameters_iter={}
	extracted_parameters_iter={}
	j=0

	while i<10:

		# Printing the iteration number
		i+=1
		print('----- Iteration ',i,' -----')

		initial_circuit_parameters=cir.get_initial_circuit_parameters()
		circuit_parameters=cir.get_circuit_parameters()
		extracted_parameters=cir.get_extracted_parameters()

		# Updating Cd
		initial_circuit_parameters['Cd']=updating_Cd(extracted_parameters['cgd2'],Cload,initial_circuit_parameters['Ld'],fo)
		
		# Updating W to improve the Qin
		Z_max=calculate_Zim_max(optimization_input_parameters['output_conditions']['s11_db'])
		Z_diff=np.abs(extracted_parameters['0_Zin_I']-extracted_parameters['2_Zin_I'])
		initial_circuit_parameters['W']=initial_circuit_parameters['W']*Z_diff/Z_max*1.2
		
		# Calculating Io from W and gm based on NF Calculation
		global gm
		gm=update_gm(extracted_parameters['nf_db'],nf,gm)
		initial_circuit_parameters['Io']=calculate_Io(gm,un,Cox,initial_circuit_parameters['W'],Lmin)

		# Running the circuit and updating the results
		cir.update_circuit(initial_circuit_parameters)
		extracted_parameters=cir.get_extracted_parameters()
		update_parameters(cir,optimization_input_parameters,i,'Ld_W_Io')
		
		# Storing the results
		initial_circuit_parameters_iter[j]=initial_circuit_parameters
		circuit_parameters_iter[j]=cir.get_circuit_parameters()
		extracted_parameters_iter[j]=extracted_parameters
		j+=1
		
		# Updating the values
		fo=optimization_input_parameters['output_conditions']['wo']/(2*np.pi)
		initial_circuit_parameters['Ls']=circuit_parameters['Ls']*50*4/(2*extracted_parameters['1_Zin_R']+extracted_parameters['0_Zin_R']+extracted_parameters['2_Zin_R'])
		initial_circuit_parameters['Lg']=initial_circuit_parameters['Lg']-(2*extracted_parameters['1_Zin_I']+extracted_parameters['0_Zin_I']+extracted_parameters['2_Zin_I'])/(4*2*np.pi*fo)
		
		# Running the circuit and updating the results
		cir.update_circuit(initial_circuit_parameters)
		update_parameters(cir,optimization_input_parameters,i,'Ls_Lg')

		# Storing the results
		initial_circuit_parameters_iter[j]=initial_circuit_parameters
		circuit_parameters_iter[j]=cir.get_circuit_parameters()
		extracted_parameters_iter[j]=cir.get_extracted_parameters()
		j+=1
	
	i=get_best_point(circuit_parameters_iter,extracted_parameters_iter,optimization_input_parameters['output_conditions'])
	
	cir.update_circuit_state(initial_circuit_parameters_iter[i],circuit_parameters_iter[i],extracted_parameters_iter[i])
	

"""
===========================================================================================================================
-------------------------------------------- Output Functions -------------------------------------------------------------
"""


#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the initial parameters by completing all the sub steps of pre optimization
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
		
	#======================================================== Step 1 =======================================================
	print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	calculate_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['auto_hc']={}
	optimization_results['auto_hc']['initial_circuit_parameters']=cir.get_initial_circuit_parameters()
	optimization_results['auto_hc']['circuit_parameters']=cir.get_circuit_parameters()
	optimization_results['auto_hc']['extracted_parameters']=cir.get_extracted_parameters()

	# Printing the values
	cf.print_initial_circuit_parameters(cir.get_initial_circuit_parameters())
	cf.print_circuit_parameters(cir.get_circuit_parameters())
	cf.print_extracted_parameters(cir.get_extracted_parameters())

	"""	
	#======================================================== Step 2 =======================================================
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	update_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['initial_circuit_parameters']=cir.get_initial_circuit_parameters()
	optimization_results['hc_update']['circuit_parameters']=cir.get_circuit_parameters()
	optimization_results['hc_update']['extracted_parameters']=cir.get_extracted_parameters()

	# Printing the values
	cf.print_initial_circuit_parameters(cir.get_initial_circuit_parameters())
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)
	"""

#===========================================================================================================================
