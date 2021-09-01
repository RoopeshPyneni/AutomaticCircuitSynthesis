#===========================================================================================================================
"""
Name: Pyneni Roopesh
Roll Number: EE18B028

Hand Calculations 1 File:
"""
#===========================================================================================================================
import numpy as np
import common_functions as cf
import spectre as sp

#===========================================================================================================================
#------------------------------------Defining the functions for simple calculations-----------------------------------------
	
#-----------------------------------------------------------------------------------------------
def update_W_Io_gm(circuit_parameters,gmRs):
	
	# Checking if the change is very high
	k=1.0/gmRs
	if k>2:
		k=2
	if k<-2:
		k=-2
	W=circuit_parameters['W']
	Io=circuit_parameters['Io']
	
	# Updating the value of W and Io
	param1=0.0 # This parameters will give the weight for whether W or Io will increase
	circuit_parameters['W']=W*(k**param1)
	circuit_parameters['Io']=Io*(k**(2-param1))
	
	return circuit_parameters

#-----------------------------------------------------------------------------------------------
def update_W_Io_gm_vdsat(circuit_parameters,gm,Rs,vdsat,vdsat_reqd):
	
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

	# Assigning the values
	L=mos_parameters['Lmin']
	cox=mos_parameters['cox']
	wo=opt_conditions['wo']
	W=circuit_parameters['W']
	
	# Calculating C2
	C2=threshold2*(2*W*L*cox/3)
	
	# Calculating Rbias
	Rbias=threshold3/(wo*C2)
	
	return C2,Rbias
	
#-----------------------------------------------------------------------------------------------
# Calculating C2 from Cgs & Cgd and Rbias
def calc_C2_updated(extracted_parameters,opt_conditions,circuit_parameters,optimization_input_parameters):
	
	threshold2=optimization_input_parameters['pre_optimization']['C2_threshold']
	threshold3=optimization_input_parameters['pre_optimization']['Rbias_threshold']

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
	

	

	
#===========================================================================================================================
#-------------------------------------------- Main Functions ---------------------------------------------------------------
	
#---------------------------------------------------------------------------------------------------------------------------
# Function to calculate the Initial Circuit Parameters	
def calculate_initial_parameters(mos_parameters,optimization_input_parameters):

	opt_conditions=optimization_input_parameters['output_conditions']
	vdsat_reqd=optimization_input_parameters['pre_optimization']['vdsat_reqd']
	
	# Creating dictionary of output parameters	
	circuit_parameters={}
	circuit_parameters_list=['Rb','Rd','Io','C1','C2','W','Rbias']
	for param_name in circuit_parameters_list:
		circuit_parameters[param_name]=0
	
	# Calculating Io min and W
	#circuit_parameters['Io'],circuit_parameters['W']=calc_Io_min(opt_conditions,mos_parameters)
	circuit_parameters['Io'],circuit_parameters['W']=calc_Io_W(opt_conditions,mos_parameters,vdsat_reqd)
	
	# Calculating Rd
	circuit_parameters['Rd']=calc_Rd(opt_conditions)
	
	# Calculating C1
	circuit_parameters['C1']=calc_C1(opt_conditions,optimization_input_parameters)
		
	# Calculating Rb
	circuit_parameters['Rb']=calc_Rb(opt_conditions,mos_parameters,circuit_parameters)	
	
	# Calculating C2 and Rbias
	circuit_parameters['C2'],circuit_parameters['Rbias']=calc_C2(mos_parameters,opt_conditions,circuit_parameters,optimization_input_parameters)
	
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(circuit_parameters,mos_parameters,opt_conditions)
	
	# Running Eldo
	extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)
	
	return circuit_parameters,dc_outputs,extracted_parameters


#---------------------------------------------------------------------------------------------------------------------------
# Function to update the Initial Circuit Parameters	
def update_initial_parameters(circuit_parameters,mos_parameters,extracted_parameters,optimization_input_parameters):

	opt_conditions=optimization_input_parameters['output_conditions']
	limit=optimization_input_parameters['pre_optimization']['Step1b_Limit']
	
	mos_parameters['vt']=extracted_parameters['vt']
	i=0
	
	while i<limit:
		print('---------------Iteration ',i+1,' --------------------')
		
		# Calculating Rb
		circuit_parameters['Rb']=calc_Rb(opt_conditions,mos_parameters,circuit_parameters)
	
		# Calculating C2 and Rbias
		circuit_parameters['C2'],circuit_parameters['Rbias']=calc_C2_updated(extracted_parameters,opt_conditions,circuit_parameters,optimization_input_parameters)
		
		# Running Eldo
		extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)
		
		# Updating the value of vt
		mos_parameters['vt']=extracted_parameters['vt']
		
		i+=1
		
	# Calculating dc outputs
	dc_outputs=calc_dc_opt(circuit_parameters,mos_parameters,opt_conditions)
	
	return circuit_parameters,dc_outputs,mos_parameters,extracted_parameters

#---------------------------------------------------------------------------------------------------------------------------
# Function to optimize gm	
def dc_optimize_gm(mos_parameters,circuit_parameters,extracted_parameters,optimization_input_parameters):
	
	# Getting the output conditions
	output_conditions=optimization_input_parameters['output_conditions']
	
	i=0
	
	limit1=optimization_input_parameters['pre_optimization']['Step2_Limit']
	 
	while i<limit1:
	
		Rs=output_conditions['Rs']
		gm=extracted_parameters['gm1']
		gmRs=gm*Rs
		
		# Checking the threshold
		threshold1=optimization_input_parameters['pre_optimization']['gmrs_threshold'] # This is the threshold for gm*Rs
		if gmRs<(1+threshold1) and gmRs>(1-threshold1):
			break
			
		print('------------------Iteration Number ',i+1,'----------------------')
	
		# Updating W and Io
		circuit_parameters=update_W_Io_gm(circuit_parameters,gmRs)
		
		# Updating Rb
		circuit_parameters['Rb']=calc_Rb(output_conditions,mos_parameters,circuit_parameters)
		
		# Running Eldo
		extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)
		
		i+=1
		
	return circuit_parameters,extracted_parameters
	

#---------------------------------------------------------------------------------------------------------------------------
# Function to optimize gm and vdsat
def dc_optimize_gm_vdsat(mos_parameters,circuit_parameters,extracted_parameters,optimization_input_parameters):
	
	# Getting the output conditions
	output_conditions=optimization_input_parameters['output_conditions']
	
	i=0
	vdsat_reqd=optimization_input_parameters['pre_optimization']['vdsat_reqd']
	
	limit1=optimization_input_parameters['pre_optimization']['Step2_Limit']
	
	while i<limit1:
	
		Rs=output_conditions['Rs']
		gm=extracted_parameters['gm1']		
		vdsat=extracted_parameters['vdsat']
		
		gmRs=gm*Rs
		
		# Checking the threshold
		threshold1=optimization_input_parameters['pre_optimization']['gmrs_threshold'] # This is the threshold for gm*Rs
		threshold2=optimization_input_parameters['pre_optimization']['vdsat_threshold'] # This is the threshold for vdsat

		if gmRs<(1+threshold1) and gmRs>(1-threshold1) and vdsat>(vdsat_reqd-threshold2) and vdsat<(vdsat_reqd+threshold2):
			break
			
		print('------------------Iteration Number ',i+1,'----------------------')
	
		# Updating W and Io
		circuit_parameters=update_W_Io_gm_vdsat(circuit_parameters,gm,Rs,vdsat,vdsat_reqd)
		
		# Updating Rb
		circuit_parameters['Rb']=calc_Rb(output_conditions,mos_parameters,circuit_parameters)
		
		# Running Eldo
		extracted_parameters=sp.write_extract(circuit_parameters,optimization_input_parameters)
		
		i+=1
		
	return circuit_parameters,extracted_parameters

#===========================================================================================================================
#-------------------------------------------- Output Functions -------------------------------------------------------------

def automatic_initial_parameters(mos_parameters,optimization_input_parameters,optimization_results):
	
	#======================================================== Step 1 =============================================================================================================

	print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Automatic Operating Point Selection 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	print('\n\n--------------------------------- Operating Point Calculations ------------------------------------')

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#--------------------Initial Point Calculations-------------------------

	# Calculating the Values of Circuit Parameters
	circuit_parameters,dc_initial_outputs,extracted_parameters=calculate_initial_parameters(mos_parameters,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['auto_hc']={}
	optimization_results['auto_hc']['circuit_parameters']=circuit_parameters.copy()
	optimization_results['auto_hc']['extracted_parameters']=extracted_parameters.copy()

	# Printing the values
	cf.print_circuit_parameters(circuit_parameters)
	cf.print_DC_outputs(dc_initial_outputs,mos_parameters)
	cf.print_extracted_outputs(extracted_parameters)

	#cf.wait_key()

	#======================================================== Step 1 b ============================================================================================================
	
	print('\n\n--------------------------------- Operating Point Updations ------------------------------------')


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#--------------------Initial Point Updating-----------------------------

	# Calculating the Values of Circuit Parameters
	circuit_parameters,dc_initial_outputs,mos_parameters,extracted_parameters=update_initial_parameters(circuit_parameters,
	mos_parameters,extracted_parameters,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['hc_update']={}
	optimization_results['hc_update']['circuit_parameters']=circuit_parameters.copy()
	optimization_results['hc_update']['extracted_parameters']=extracted_parameters.copy()

	# Printing the values
	cf.print_circuit_parameters(circuit_parameters)
	cf.print_DC_outputs(dc_initial_outputs,mos_parameters)
	cf.print_extracted_outputs(extracted_parameters)

	#cf.wait_key()


	#======================================================== Step 2 =============================================================================================================

	print('\n\n--------------------------------- gm and vdsat Updation ------------------------------------')
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#--------------------gm and vdsat updating------------------------------
	circuit_parameters,extracted_parameters=dc_optimize_gm_vdsat(mos_parameters,circuit_parameters,extracted_parameters,optimization_input_parameters)

	# Storing the Circuit and Extracted Parameters
	optimization_results['gmvd_update']={}
	optimization_results['gmvd_update']['circuit_parameters']=circuit_parameters.copy()
	optimization_results['gmvd_update']['extracted_parameters']=extracted_parameters.copy()

	# Printing the values
	cf.print_circuit_parameters(circuit_parameters)
	cf.print_extracted_outputs(extracted_parameters)

	#cf.wait_key()

	return circuit_parameters,extracted_parameters

#===========================================================================================================================


		
	


