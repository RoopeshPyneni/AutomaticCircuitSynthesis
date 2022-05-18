#===========================================================================================================================
"""
Name				: Roopesh Pyneni
Roll Number			: EE18B028
File Description 	: This file will perform hand calculations using method 1
					  This method involves choosing the points such that gm=20mS and Rb is high
"""

#===========================================================================================================================
import numpy as np
import os
import common_functions as cf # type: ignore

"""
===========================================================================================================================
--------------------------------------------- Defining the functions for simple calculations ------------------------------
"""

#-----------------------------------------------------------------------------------------------
# This function will update the value of W and Io so that gm improves
def update_W_Io_gm(circuit_parameters,gmRs):
	
	# Calculating the value by which we should change W and Io
	k=1.0/gmRs
	
	# Checking if the change is very high
	# If it is high, we will limit the change to 2
	if k>2:
		k=2
	if k<-2:
		k=-2
	W=circuit_parameters['W']
	Io=circuit_parameters['Io']
	
	# This parameters will give the weight for whether W or Io will increase
	# If this is 0, we will not update W  and will only update Io
	# If this is 2, we will not update IO and will only update W 
	# Choose a value between 0 and 2
	param1=0.0 

	# Updating the value of W and Io
	circuit_parameters['W']=W*(k**param1)
	circuit_parameters['Io']=Io*(k**(2-param1))
	
	return circuit_parameters

#-----------------------------------------------------------------------------------------------
# This function will update the value of W and Io so that gm and vdsat improves
def update_W_Io_gm_vdsat(circuit_parameters,gm,Rs,vdsat,vdsat_reqd):
	
	# Calculating the value by which we should change W and Io
	k1=1.0/(gm*Rs)
	k2=vdsat_reqd/vdsat
	
	# Checking if the change is very high
	if k1>2:
		k1=2
	
	if k2>2:
		k2=2
	if k2<0:
		k2=2
		
	# Updating the value of W and Io
	circuit_parameters['W']=circuit_parameters['W']*k1/k2
	circuit_parameters['Io']=circuit_parameters['Io']*k1*k2
	
	return circuit_parameters

#-----------------------------------------------------------------------------------------------
# Calculating Io min and W max by considering the s11 parameter
def calc_Io_min(opt_conditions,mos_parameters):

	# Assigning the values
	s11_db=opt_conditions['s11_db']
	wo=opt_conditions['wo']
	Rs=opt_conditions['Rs']
	Lmin=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	un=mos_parameters['un']
	
	# Calculating Y
	s11_normal=cf.db_to_normal(s11_db)
	Y=(2/Rs)*np.sqrt(s11_normal)
	
	# Calculating Wmax and gm
	wmax=3*Y/(2*wo*Lmin*cox)
	gm=1/Rs

	# Calculating Io,min	
	Io_min=(gm*gm*Lmin)/(2*un*cox*wmax)
	
	return Io_min,wmax
	
#-----------------------------------------------------------------------------------------------
# Calculating Io and W max by considering vd sat
def calc_Io_W(opt_conditions,mos_parameters,vdsat_reqd):

	# Assigning the values
	Rs=opt_conditions['Rs']
	Lmin=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	un=mos_parameters['un']
	
	# Calculating W
	W=Lmin/(Rs*un*cox*vdsat_reqd)
	
	# Calculating Io	
	Io=vdsat_reqd/(2*Rs)
	
	return Io,W

#-----------------------------------------------------------------------------------------------	
# Calculating Rd by assuming gain=Rd/2*Rs
def calc_Rd(opt_conditions):

	# Calculating the required gain
	gain=np.sqrt(cf.db_to_normal(opt_conditions['gain_db']))

	# Assigning the values
	Rs=opt_conditions['Rs']
	nf_db=opt_conditions['nf_db']
	
	# Calculating f from nf
	f=cf.db_to_normal(nf_db)
	
	# Assigning value for gamma
	gamma=2
	
	# Calculating Rd
	Rd1=2*Rs*gain
	Rd2=5*Rs/(f-1-gamma)
	Rd=np.maximum(Rd1,Rd2)
	
	return Rd

#-----------------------------------------------------------------------------------------------
# Calculating C1
def calc_C1(opt_conditions,optimization_input_parameters):
	
	threshold1=optimization_input_parameters['pre_optimization']['C1_threshold']

	# Assigning the values
	wo=opt_conditions['wo']
	Rs=opt_conditions['Rs']
	
	# Calculating C1
	C1=threshold1/(wo*Rs)
	
	return C1

#-----------------------------------------------------------------------------------------------
# Calculating Rb
def calc_Rb(opt_conditions,mos_parameters,circuit_parameters):

	# Assigning the values
	vdd=mos_parameters['vdd']
	vt=mos_parameters['vt']
	Rs=opt_conditions['Rs']
	delta_v=opt_conditions['delta_v']
	iip3_dbm=opt_conditions['iip3_dbm']
	Io=circuit_parameters['Io']
	Rd=circuit_parameters['Rd']
	
	# Calculating gain
	gain=Rd/(2*Rs)

	# Calculating p1,db
	p1_db=iip3_dbm-39.6
	
	# Calculating the input swing
	vi_swing=np.sqrt(8*Rs*cf.db_to_normal(p1_db))
	
	# Calculating the output swing
	vo_swing_min=vi_swing*gain
	
	# Calculating vd sat
	vdsat=2*Io*Rs	# Here, gm=1/Rs
	
	# Calculating Rb1
	Rb1=(vdd-Io*Rd-vdsat-vo_swing_min)/Io
	
	# Calculating Rb2
	Rb2=(vdd-delta_v-vdsat-vt)/Io
	
	# Calculating Rb
	Rb=np.minimum(Rb1,Rb2)
	
	if Rb<10:
		Rb=10
	
	return Rb

#-----------------------------------------------------------------------------------------------
# Calculating C2 from W and Rbias 
def calc_C2(mos_parameters,opt_conditions,circuit_parameters,optimization_input_parameters):
	
	threshold2=optimization_input_parameters['pre_optimization']['C2_threshold']
	threshold3=optimization_input_parameters['pre_optimization']['Rbias_threshold']
	Rbias_min=optimization_input_parameters['pre_optimization']['Rbias_minimum']

	# Assigning the values
	L=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	wo=opt_conditions['wo']
	W=circuit_parameters['W']
	
	# Calculating C2
	C2=threshold2*(2*W*L*cox/3)
	
	# Calculating Rbias
	Rbias=threshold3/(wo*C2)
	Rbias=max(Rbias,Rbias_min)
	
	return C2,Rbias
	
#-----------------------------------------------------------------------------------------------
# Calculating C2 from Cgs & Cgd and Rbias
def calc_C2_updated(extracted_parameters,opt_conditions,circuit_parameters,optimization_input_parameters):
	
	threshold2=optimization_input_parameters['pre_optimization']['C2_threshold']
	threshold3=optimization_input_parameters['pre_optimization']['Rbias_threshold']
	Rbias_min=optimization_input_parameters['pre_optimization']['Rbias_minimum']

	# Assigning the values
	cgs=extracted_parameters['cgs1']
	cgd=extracted_parameters['cgd1']
	gain=circuit_parameters['Rd']/(2*opt_conditions['Rs'])
	wo=opt_conditions['wo']
	
	# Calculating C2
	C2a=threshold2*cgs
	C2b=threshold2*cgd*gain
	C2=np.maximum(C2a,C2b)
	
	# Calculating Rbias
	Rbias=threshold3/(wo*C2)
	Rbias=max(Rbias,Rbias_min)
	
	return C2,Rbias
	
#-----------------------------------------------------------------------------------------------
# Calculating DC Output Values from Initial Circuit Parameters
def calc_dc_opt(circuit_parameters,mos_parameters,opt_conditions):
	vs=circuit_parameters['Io']*circuit_parameters['Rb']
	vgs=mos_parameters['vt']+ 2*opt_conditions['Rs']*circuit_parameters['Io']
	vg=vgs+vs
	vd=mos_parameters['vdd']-circuit_parameters['Io']*circuit_parameters['Rd']
	
	dc_outputs={'vg':vg,'vd':vd,'vs':vs}
	
	return dc_outputs
	

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

	mos_parameters=cir.mos_parameters
	opt_conditions=optimization_input_parameters['output_conditions']
	vdsat_reqd=optimization_input_parameters['pre_optimization']['vdsat_reqd']
	
	# Creating dictionary of output parameters	
	initial_circuit_parameters={}
	circuit_parameters_list=['Rb','Rd','Io','C1','C2','W','Rbias']
	for param_name in circuit_parameters_list:
		initial_circuit_parameters[param_name]=0
	
	# Calculating Io min and W
	#initial_circuit_parameters['Io'],initial_circuit_parameters['W']=calc_Io_min(opt_conditions,mos_parameters)
	initial_circuit_parameters['Io'],initial_circuit_parameters['W']=calc_Io_W(opt_conditions,mos_parameters,vdsat_reqd)
	print(initial_circuit_parameters)
	
	# Calculating Rd
	initial_circuit_parameters['Rd']=calc_Rd(opt_conditions)
	
	# Calculating C1
	initial_circuit_parameters['C1']=calc_C1(opt_conditions,optimization_input_parameters)
		
	# Calculating Rb
	initial_circuit_parameters['Rb']=calc_Rb(opt_conditions,mos_parameters,initial_circuit_parameters)	
	
	# Calculating C2 and Rbias
	initial_circuit_parameters['C2'],initial_circuit_parameters['Rbias']=calc_C2(mos_parameters,opt_conditions,initial_circuit_parameters,optimization_input_parameters)
	
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(initial_circuit_parameters,mos_parameters,opt_conditions)
	
	# Running Eldo
	cir.update_circuit(initial_circuit_parameters)
	
	return dc_outputs


#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters	after calculating the new value of vt
def update_initial_parameters(cir,optimization_input_parameters):

	opt_conditions=optimization_input_parameters['output_conditions']
	limit=optimization_input_parameters['pre_optimization']['Step1b_Limit']
	
	cir.mos_parameters['vt']=cir.extracted_parameters['vt']
	i=0
	
	while i<limit:
		print('---------------Iteration ',i+1,' --------------------')
		
		# Calculating Rb
		cir.initial_circuit_parameters['Rb']=calc_Rb(opt_conditions,cir.mos_parameters,cir.initial_circuit_parameters)
	
		# Calculating C2 and Rbias
		cir.initial_circuit_parameters['C2'],cir.initial_circuit_parameters['Rbias']=calc_C2_updated(cir.extracted_parameters,opt_conditions,cir.initial_circuit_parameters,optimization_input_parameters)
		
		# Running Eldo
		cir.run_circuit()


		
		# Updating the value of vt
		cir.mos_parameters['vt']=cir.extracted_parameters['vt']
		
		i+=1
		
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(cir.initial_circuit_parameters,cir.mos_parameters,opt_conditions)
	
	return dc_outputs

#---------------------------------------------------------------------------------------------------------------------------
# Function to change circuit parameters to get better gm
def dc_optimize_gm(cir,optimization_input_parameters):
	
	# Getting the output conditions
	output_conditions=optimization_input_parameters['output_conditions']
	
	i=0
	
	limit1=optimization_input_parameters['pre_optimization']['Step2_Limit']
	 
	while i<limit1:
	
		Rs=output_conditions['Rs']
		gm=cir.extracted_parameters['gm1']
		gmRs=gm*Rs
		
		# Checking the threshold
		threshold1=optimization_input_parameters['pre_optimization']['gmrs_threshold'] # This is the threshold for gm*Rs
		if gmRs<(1+threshold1) and gmRs>(1-threshold1):
			break
			
		print('------------------Iteration Number ',i+1,'----------------------')
	
		# Updating W and Io
		cir.initial_circuit_parameters=update_W_Io_gm(cir.initial_circuit_parameters,gmRs)
		
		# Updating Rb
		cir.initial_circuit_parameters['Rb']=calc_Rb(output_conditions,cir.mos_parameters,cir.initial_circuit_parameters)
		
		# Running Eldo
		cir.run_circuit()
		
		i+=1
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to optimize gm and vdsat
def dc_optimize_gm_vdsat(cir,optimization_input_parameters):
	
	# Getting the output conditions
	output_conditions=optimization_input_parameters['output_conditions']
	
	i=0
	vdsat_reqd=optimization_input_parameters['pre_optimization']['vdsat_reqd']
	
	limit1=optimization_input_parameters['pre_optimization']['Step2_Limit']
	
	while i<limit1:
	
		Rs=output_conditions['Rs']
		gm=cir.extracted_parameters['gm1']		
		vdsat=cir.extracted_parameters['vdsat']
		
		gmRs=gm*Rs
		
		# Checking the threshold
		threshold1=optimization_input_parameters['pre_optimization']['gmrs_threshold'] # This is the threshold for gm*Rs
		threshold2=optimization_input_parameters['pre_optimization']['vdsat_threshold'] # This is the threshold for vdsat

		if gmRs<(1+threshold1) and gmRs>(1-threshold1) and vdsat>(vdsat_reqd-threshold2) and vdsat<(vdsat_reqd+threshold2):
			break
			
		print('------------------Iteration Number ',i+1,'----------------------')
	
		# Updating W and Io
		cir.initial_circuit_parameters=update_W_Io_gm_vdsat(cir.initial_circuit_parameters,gm,Rs,vdsat,vdsat_reqd)
		
		# Updating Rb
		cir.initial_circuit_parameters['Rb']=calc_Rb(output_conditions,cir.mos_parameters,cir.initial_circuit_parameters)
		
		# Running Eldo
		cir.run_circuit()
		
		i+=1


"""
===========================================================================================================================
-------------------------------------------- Output Functions -------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the initial parameters by completing all the sub steps of pre optimization
def automatic_initial_parameters(cir,optimization_input_parameters,optimization_results):
		
	print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


	#======================================================== Step 1 ========================================================
	print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	dc_initial_outputs=calculate_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['auto_hc']={}
	optimization_results['auto_hc']['initial_circuit_parameters']=cir.initial_circuit_parameters.copy()
	optimization_results['auto_hc']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['auto_hc']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cf.print_initial_circuit_parameters(cir.initial_circuit_parameters)
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)

	

	#======================================================== Step 1 b ======================================================
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')

	# Calculating the Values of Circuit Parameters
	dc_initial_outputs=update_initial_parameters(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['initial_circuit_parameters']=cir.initial_circuit_parameters.copy()
	optimization_results['hc_update']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['hc_update']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cf.print_initial_circuit_parameters(cir.initial_circuit_parameters)
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)



	#======================================================== Step 2 ========================================================
	print('\n\n--------------------------------- gm and vdsat Updation ------------------------------------')
	
	# Calculating the Values of Circuit Parameters
	dc_optimize_gm_vdsat(cir,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['gmvd_update']={}
	optimization_results['gmvd_update']['initial_circuit_parameters']=cir.initial_circuit_parameters.copy()
	optimization_results['gmvd_update']['circuit_parameters']=cir.circuit_parameters.copy()
	optimization_results['gmvd_update']['extracted_parameters']=cir.extracted_parameters.copy()

	# Printing the values
	cf.print_initial_circuit_parameters(cir.initial_circuit_parameters)
	cf.print_circuit_parameters(cir.circuit_parameters)
	cf.print_extracted_parameters(cir.extracted_parameters)

#===========================================================================================================================
