#!/usr/bin/env python -i


import csv, sys, os, array, warnings
import input_sub as input_sub
import numpy as np
from astropy import constants as c
import pyPLUTO as pp
from matplotlib import pyplot as plt
from astropy.io import ascii

gamma=5./3.


PSD_density=ascii.read("PSD_8_density.csv")
PSD_vr=ascii.read("PSD_8_v_r.csv")
PSD_rho_vr=ascii.read("PSD_8_rho_vr.csv")
PSD_dmdot=ascii.read("PSD_8_dmdot.csv")



#First get scaling factorz from the definitions file

inp=open('definitions.h','r')
for line in inp.readlines():
	data=line.split()
	if len(data)>1:
		if data[1]=='UNIT_DENSITY':
			UNIT_DENSITY=float(data[2])
		elif data[1]=='UNIT_LENGTH':
			UNIT_LENGTH=float(data[2])
		elif data[1]=='UNIT_VELOCITY':
			UNIT_VELOCITY=float(data[2])

#Compute deived scaling factors

UNIT_MASS=(UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
UNIT_ACCELERATION=(UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
UNIT_FORCE=(UNIT_MASS*UNIT_ACCELERATION)
UNIT_TIME=(UNIT_LENGTH/UNIT_VELOCITY)
UNIT_PRESSURE=(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)

#Compute the numebr that transforms from pressure to temperature


KELVIN=UNIT_VELOCITY*UNIT_VELOCITY*c.m_p.cgs/c.k_B.cgs


inp.close()


#open the actual data file

fname=int(sys.argv[1])
D=pp.pload(fname)


#Get the pressure and density

density=np.transpose(D.rho)*UNIT_DENSITY

try:
    pressure=np.transpose(D.prs)*UNIT_PRESSURE
    energy=pressure/(gamma-1.)
    E_int=pressure/(2./3.)
    
    
except:
    print ("No pressure info - must be isothermal")

#Compute the internal energy from the pressure


#Compute/get number densities

#nd=density/(1.43*c.m_p.value)
try:
	ne=np.transpose(D.ne)
	nh=np.transpose(D.nh)
except:
	print("No ne or nh fields, using 1.43 as scaling to nh")
	nh=density/(1.43*c.m_p.cgs).value
	ne=nh*1.21

#Get the velocities

v_r=np.transpose(D.vx1)*UNIT_VELOCITY
v_t=np.transpose(D.vx2)*UNIT_VELOCITY
v_p=np.transpose(D.vx3)*UNIT_VELOCITY

E_kin=density*(v_r**2.+v_t**2.+v_p**2.)

#And compute the speed

v=np.sqrt(v_r**2+v_t**2)

#Get the cooling rates (if here)

try:
	line_c=np.transpose(D.line_c)
	xray_h=np.transpose(D.xray_h)
	comp_c=np.transpose(D.comp_c)
	comp_h=np.transpose(D.comp_h)
	brem_c=np.transpose(D.brem_c)
except:
	print("No old style cooling rate info ")
	
try:
	line_c=np.transpose(D.lc)
	xray_h=np.transpose(D.xh)
	comp_c=np.transpose(D.cc)
	comp_h=np.transpose(D.ch)
	brem_c=np.transpose(D.bc)
	print("Read cooling rates as line_c,xray_h,comp_c,comp_h,brem_c")
    
except:
	print("No new style cooling rate info")	
	

try:
	line_c_pre=np.transpose(D.line_c_pre)
	xray_h_pre=np.transpose(D.xray_h_pre)
	comp_c_pre=np.transpose(D.comp_c_pre)
	comp_h_pre=np.transpose(D.comp_h_pre)
	brem_c_pre=np.transpose(D.brem_c_pre)
except:
	print("No old style cooling rate prefactors")
	
	
try:
	line_c_pre=np.transpose(D.lc_pre)
	xray_h_pre=np.transpose(D.xh_pre)
	comp_c_pre=np.transpose(D.cc_pre)
	comp_h_pre=np.transpose(D.ch_pre)
	brem_c_pre=np.transpose(D.bc_pre)
	print("Read cooling rate prefactors as line_c_pre,xray_h_pre,comp_c_pre,comp_h_pre,brem_c_pre")
except:
	print("No new style cooling rate prefactors")	
	
	
try:
	g_r=np.transpose(D.gr)*UNIT_ACCELERATION
	g_t=np.transpose(D.gt)*UNIT_ACCELERATION
	g_p=np.transpose(D.gp)*UNIT_ACCELERATION
	print("Read rad driving accelerations as g_r,g_t,g_p")
except:
	print("No Radiation acceleration")
	
try:
	gx_pre=np.transpose(D.gx_pre)
	gy_pre=np.transpose(D.gy_pre)
	gz_pre=np.transpose(D.gz_pre)
	print("Read rad driving prefactors as gx_pre,gy_pre,gz_pre")
except:
	print("No Radiation acceleration prefactors")	
	
	
try:
	tau_es=np.transpose(D.t_es)
except:
	print("No electron sacattering optcial depth")

#Get optcially thin ionization parameter if here

try:
	xi=np.transpose(D.XI)
except:
	print("No ionization parameter")
	
#Get temperature - if present or calculate it
	
try:
	temperature=np.transpose(D.T)
except:
    try:
	    print("No temperature data - computing")
	    temperature=pressure/UNIT_PRESSURE*KELVIN*0.6/(density/UNIT_DENSITY)
    except:
        print("no pressure info")
        
    
try:
    ciso=np.sqrt(5./3.*c.k_B.cgs.value*temperature/c.m_p.cgs.value)    
except:
    print("no pressure info")    
    

#Get the geometric quantities

r=D.x1*UNIT_LENGTH
theta=D.x2

#See if we have thermal dt

try:
	dt=np.transpose(D.dt)
except:
	print("No thermal timestep info")

plt.plot(np.degrees(theta),v_r[::,-1]/1000./100.,label="PY PLUTO")
plt.plot(PSD_vr["col1"],PSD_vr["col2"],label="PSD")
plt.xlabel("theta")
plt.ylabel("v_r (km/s)")
plt.legend()

plt.ylim(0,20000)
plt.tight_layout()
plt.savefig("zzzz_"+sys.argv[1]+"_v_r_zoom.png")
plt.close()

plt.plot(np.degrees(theta),v_r[::,-1]/1000./100.,label="PY PLUTO")
plt.plot(PSD_vr["col1"],PSD_vr["col2"],label="PSD")
plt.xlabel("theta")
plt.ylabel("v_r (km/s)")
plt.legend()
plt.tight_layout()
plt.savefig("zzzz_"+sys.argv[1]+"_v_r.png")
plt.close()


plt.semilogy(np.degrees(theta),density[::,-1],label="PY PLUTO")
plt.plot(PSD_density["col1"],PSD_density["col2"],label="PSD")
plt.xlabel("theta")
plt.ylabel("density (g/cm^3)")
plt.legend()
plt.tight_layout()
plt.savefig("zzzz_"+sys.argv[1]+"_rho.png")
plt.close()

plt.semilogy(np.degrees(theta),v_r[::,-1]*density[::,-1],label="PY PLUTO")
plt.plot(PSD_rho_vr["col1"],PSD_rho_vr["col2"],label="PSD")
plt.xlabel("theta")
plt.ylabel("density*v_r (g/cm^2/s)")
plt.legend()
plt.tight_layout()
plt.savefig("zzzz_"+sys.argv[1]+"_rhovr.png")
plt.close()





theta_ratio=(theta[2]-theta[1])/(theta[1]-theta[0])
dtheta=[]
dtheta.append((theta[1]-theta[0])/(0.5*(1.0+theta_ratio)))
theta_min=theta[0]-0.5*dtheta[-1]
theta_max=theta_min+dtheta[-1]*(1.0-theta_ratio**len(theta))/(1.0-theta_ratio)
for i in range(len(theta)-1):
	dtheta.append(theta_ratio*dtheta[-1])
    
    
dmdot=[]
cum_dmdot=[]
for i in range(len(theta)):
    dmdot.append(density[::,-1][i]*v_r[::,-1][i]*np.sin(theta[i])*dtheta[i])   
    if i==0:
        cum_dmdot.append(dmdot[-1])
    else:
        cum_dmdot.append(dmdot[-1]+cum_dmdot[-1])
    
dmdot=np.array(dmdot)*4.*np.pi*r[-1]*r[-1]     
cum_dmdot=np.array(cum_dmdot)*4.*np.pi*r[-1]*r[-1]     


plt.semilogy(np.degrees(theta),cum_dmdot,label="PY PLUTO")
plt.semilogy(PSD_dmdot["col1"],PSD_dmdot["col2"],label="PSD")
plt.ylim(1e6,1e16)
plt.xlabel("theta")
plt.ylabel("dmdot_theta  (g/s)")
plt.legend()
plt.tight_layout()
plt.savefig("zzzz_"+sys.argv[1]+"_dmdot.png")
plt.close()