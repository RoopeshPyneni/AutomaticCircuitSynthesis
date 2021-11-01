import numpy as np

#-----------------------------------------------------------------------------------------------
# Converting from dB to normal
# Outputs : normal values
def db_to_normal(val_db):
    return 10**(val_db/10)

#-----------------------------------------------------------------------------------------------
# Calculating the resistance from the inductance values
# Outputs : resistance
def calculate_resistance_inductor(ind,fo):
    q=15
    res=q*2*np.pi*fo*ind
    return res

#-----------------------------------------------------------------------------------------------
# Calculating Ld and Rd
# Outputs : Ld, Rd
def calculate_Ld_Rd(Cload,fo):
    wo=2*np.pi*fo
    Ld=1/(wo*wo*Cload)
    Rd=calculate_resistance_inductor(Ld,fo)
    return Ld,Rd

#-----------------------------------------------------------------------------------------------
# Calculating gm
# Outputs : gm
def calculate_gm(gain,Rd):
    return gain/(2*Rd*2.5)

#-----------------------------------------------------------------------------------------------
# Calculating Cgs
# Outputs : cgs
def calculate_cgs(Rs,Rd,F,fo,gm):
    wo=2*np.pi*fo
    wt=wo*np.sqrt((2*gm*Rs+4*Rs/Rd)/(F-1))
    cgs=gm/wt
    return cgs

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
# Calculating Ls and Rs
# Outputs : Ls,Rs
def calculate_Ls_Rls(rsource,cgs,gm,fo):
    Ls=rsource*cgs/gm
    Rls=calculate_resistance_inductor(Ls,fo)
    return Ls,Rls

#-----------------------------------------------------------------------------------------------
# Calculating Lg and Rg
# Outputs : Lg,Rg
def calculate_Lg_Rg(Ls,cgs,fo):
    w=2*np.pi*fo
    Lg=1/(w*w*cgs)-Ls
    Rg=calculate_resistance_inductor(Lg,fo)
    return Lg,Rg


def calculate_initial_parameters(cir,output_conditions):
    
	# Getting the output conditions
    Cload=output_conditions['Cload']
    fo=output_conditions['fo']
    Rs=output_conditions['Rs']
    gain=db_to_normal(output_conditions['gain_db'])
    F=db_to_normal(output_conditions['nf_db'])
    Lmin=cir.mos_parameters['Lmin']
    Cox=cir.mos_parameters['Cox']
    un=cir.mos_parameters['un']

    # Calculating the circuit parameters
    circuit_parameters={}
    circuit_parameters['Ld'],circuit_parameters['Rd']=calculate_Ld_Rd(Cload,fo)
    gm=calculate_gm(gain,circuit_parameters['Rd'])
    cgs=calculate_cgs(Rs,circuit_parameters['Rd'],F,fo,gm)
    circuit_parameters['W']=calculate_W(cgs,Lmin,Cox)
    circuit_parameters['Io']=calculate_Io(gm,un,Cox,circuit_parameters['W'],Lmin)
    circuit_parameters['Ls'],circuit_parameters['Rls']=calculate_Ls_Rls(Rs,cgs,gm,fo)
    circuit_parameters['Lg'],circuit_parameters['Rg']=calculate_Lg_Rg(circuit_parameters['Ls'],cgs,fo)
    circuit_parameters['Rb']=5000
    circuit_parameters['Cs']=100/(2*np.pi*50*fo)

	# Running the circuit
    cir.update_circuit(circuit_parameters)