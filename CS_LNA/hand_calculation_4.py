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


"""
===========================================================================================================================
------------------------------------Defining the functions for simple calculations-----------------------------------------
"""


#-----------------------------------------------------------------------------------------------
# Converting from dB to normal
# Outputs : normal values
def db_to_normal(val_db):
    return 10**(val_db/10)

#-----------------------------------------------------------------------------------------------
# Calculating Ld
# Outputs : Ld
def calculate_Ld(Cload,Cd,fo):
    wo=2*np.pi*fo
    Ld=1/(wo*wo*(Cload+Cd))
    return Ld

#-----------------------------------------------------------------------------------------------
# Calculating Qin
# Outputs : Qin
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
# Outputs : cgs
def calculate_cgs(Rs,fo,Qin):
    wo=2*np.pi*fo
    cgs=1/(2*Rs*Qin*wo)
    return cgs

#-----------------------------------------------------------------------------------------------
# Calculating gm
# Outputs : gm
def calculate_gm(Ld,fo,cgs,nf):
	wo=2*np.pi*fo
	Rd=wo*Ld*15
	A=200*((wo*cgs)**2)/Rd
	B=100*((wo*cgs)**2)
	C=1-10**(0.1*nf)
	return 2*A/(np.sqrt(B*B-4*A*C)-B)


#-----------------------------------------------------------------------------------------------
# Calculating W
# Outputs : W
def calculate_W(cgs,Lmin,Cox):
    return 3*cgs/(2*Lmin*Cox)

#-----------------------------------------------------------------------------------------------
# Calculating Io
# Outputs : Io
def calculate_Io(gm,un,cox,W,Lmin):
    return (gm*gm)/(2*un*cox*W/Lmin)

#-----------------------------------------------------------------------------------------------
# Calculating Ls
# Outputs : Ls
def calculate_Ls(rsource,cgs,gm,fo):
    return rsource*cgs/gm
    
#-----------------------------------------------------------------------------------------------
# Calculating Lg
# Outputs : Lg
def calculate_Lg(Ls,cgs,fo):
    w=2*np.pi*fo
    Lg=1/(w*w*cgs)-Ls
    return Lg

#-----------------------------------------------------------------------------------------------
# Calculating Zim_max
# Outputs : Zim_max
def calculate_Zim_max(s11):
	s11_mag_sq=10**(0.1*s11)
	H_by_Zim_max_sq=(1/s11_mag_sq)-1
	Z_im_max=100*np.sqrt(1/H_by_Zim_max_sq)
	return Z_im_max

#-----------------------------------------------------------------------------------------------
# Calculating R1 and R2
def calculate_R1_R2(gm,Io,vth,vdd,fo,Cg):
	Reff=10/(2*np.pi*fo*Cg)
	vg=0.9*vdd
	R1=Reff*vdd/vg
	R2=R1*vg/(vdd-vg)
	return R1,R2

#-----------------------------------------------------------------------------------------------
# Updating R1 and R2
def updating_R1_R2(vdsat1,vdsat2,vth,vdd):
	vg=vdsat1+vdsat2+vth
	R2=1e4*vg/vdd
	R1=1e4-R2
	return R1,R2

#-----------------------------------------------------------------------------------------------
# Updating Ld
def updating_Ld(Cgd,Cload,Cd,fo):
	wo=2*np.pi*fo
	Ld=1/(wo*wo*(Cload+Cd+Cgd))
	return Ld

#-----------------------------------------------------------------------------------------------
# Updating gm
def update_gm(extracted_nf,target_nf,gm):

	# Getting the value of f
	extracted_f=cf.db_to_normal(extracted_nf)
	target_f=cf.db_to_normal(target_nf)
	k=(extracted_f-1)/(target_f-1)

	if extracted_nf<(target_nf):
		return gm*np.sqrt(k)
	else:
		return gm*k
	

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
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/circuit_parameters.csv'

	f=open(filename,'w')
	f.write('Iter_Number,Iter_Type')
	for param_name in cir.circuit_parameters:
		f.write(','+param_name)
	f.write('\n')
	f.close()

	# Storing results for extracted parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/extracted_parameters.csv'

	f=open(filename,'w')
	f.write('Iter_Number,Iter_Type')
	for param_name in cir.extracted_parameters:
		f.write(','+param_name)
	f.write('\n')
	f.close()

#---------------------------------------------------------------------------------------------------------------------------
# Writing the values of circuit parameters and extracted parameters from each iteration to a csv file
def update_parameters(cir,optimization_input_parameters,iter_no,iter_type):
	
	# Storing results for circuit parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/circuit_parameters.csv'
	
	f=open(filename,'a')
	f.write(str(iter_no)+','+str(iter_type))
	for param_name in cir.circuit_parameters:
		f.write(','+str(cir.circuit_parameters[param_name]))
	f.write('\n')
	f.close()

	# Storing results for extracted parameters
	filename=optimization_input_parameters['filename']['output']+'/Pre_Optimization/HC_Update/Results/extracted_parameters.csv'
	
	f=open(filename,'a')
	f.write(str(iter_no)+','+str(iter_type))
	for param_name in cir.extracted_parameters:
		f.write(','+str(cir.extracted_parameters[param_name]))
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
	circuit_parameters={}
	circuit_parameters['Cd']=Cload
	circuit_parameters['Ld']=calculate_Ld(Cload,circuit_parameters['Cd'],fo)
	Qin=calculate_Qin(s11,fo,f_range)
	cgs=calculate_cgs(Rs,fo,Qin)
	circuit_parameters['W']=calculate_W(cgs,Lmin,Cox)
	global gm
	gm=calculate_gm(circuit_parameters['Ld'],fo,cgs,nf)
	circuit_parameters['Cg']=100*cgs
	circuit_parameters['Io']=calculate_Io(gm,un,Cox,circuit_parameters['W'],Lmin)
	circuit_parameters['R1'],circuit_parameters['R2']=calculate_R1_R2(gm,circuit_parameters['Io'],vth,vdd,fo,circuit_parameters['Cg'])
	circuit_parameters['Ls']=calculate_Ls(Rs,cgs,gm,fo)
	circuit_parameters['Lg']=calculate_Lg(circuit_parameters['Ls'],cgs,fo)
	circuit_parameters['Rb']=5000
	circuit_parameters['Cs']=100/(2*np.pi*50*fo)

	# Running the circuit
	cir.update_circuit(circuit_parameters)


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

	while i<10: #and cir.extracted_parameters['s11_db']>optimization_input_parameters['output_conditions']['s11_db']:

		# Printing the iteration number
		i+=1
		print('----- Iteration ',i,' -----')

		# Updating Ld
		cir.circuit_parameters['Ld']=updating_Ld(cir.extracted_parameters['cgd2'],Cload,cir.circuit_parameters['Cd'],fo)

		# Updating W to improve the Qin
		Z_max=calculate_Zim_max(optimization_input_parameters['output_conditions']['s11_db'])
		Z_diff=0.5*np.abs(cir.extracted_parameters['0_Zin_I']-cir.extracted_parameters['2_Zin_I'])
		cir.circuit_parameters['W']=cir.circuit_parameters['W']*Z_diff/Z_max
		
		# Calculating Io from W and gm based on NF Calculation
		global gm
		gm=update_gm(cir.extracted_parameters['nf_db'],nf,gm)
		cir.circuit_parameters['Io']=calculate_Io(gm,un,Cox,cir.circuit_parameters['W'],Lmin)

		# Running the circuit and updating the results
		cir.run_circuit()
		update_parameters(cir,optimization_input_parameters,i,'Ld_W_Io')
		
		# Updating the values
		fo=optimization_input_parameters['output_conditions']['wo']/(2*np.pi)
		cir.circuit_parameters['Ls']=cir.circuit_parameters['Ls']*50*4/(2*cir.extracted_parameters['1_Zin_R']+cir.extracted_parameters['0_Zin_R']+cir.extracted_parameters['2_Zin_R'])
		cir.circuit_parameters['Lg']=cir.circuit_parameters['Lg']-(2*cir.extracted_parameters['1_Zin_I']+cir.extracted_parameters['0_Zin_I']+cir.extracted_parameters['2_Zin_I'])/(4*2*np.pi*fo)
		
		# Running the circuit and updating the results
		cir.run_circuit()
		update_parameters(cir,optimization_input_parameters,i,'Ls_Lg')
	

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
	optimization_results['auto_hc']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['auto_hc']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)

	

	#======================================================== Step 2 =======================================================
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	update_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['hc_update']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)

#===========================================================================================================================
