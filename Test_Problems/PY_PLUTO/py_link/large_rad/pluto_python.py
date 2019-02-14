#!/usr/bin/env python

import subprocess
import glob
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO as pp
import numpy as np
import pluto_python_sub as pps



nproc_py=4
nproc_pl=4

rad_force=1



UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=pps.get_units()


#Set the parameters of the run here
data={}
#First set the temperature of the radiation and the central source mass
data["T_x"]=5.6e7
data["CENT_MASS"]=7.0*c.M_sun.cgs.value
data["MU"]=0.6

#We can now compute R_IC - helpful for scalings

R_IC=(c.G.cgs*data["CENT_MASS"]*data["MU"]*c.m_p.cgs/c.k_B.cgs/(data["T_x"]/4.0)).value

#The next three values define the midplane density 

data["RHO_ALPHA"]=2.0    #The drop off with radius
data["R_0"]=R_IC         #THe radius we are using to set the density (usually R_IC)
data["RHO_0"]=16.e-12   #The density at that radius

#Now lets set up the central source - we want everything to be consistent!

efficiency=0.083        #The efficiency of conversion of mass to lumonisity at the central source - used to set the disk mdot
data["L_x"]=3.3e37       #The luminosity from 13.6eV to infinity
data["BREM_ALPHA"]=0.0   #The power law for the bremstrahlung specrrum - should stay at zero unless there is a very goo reason
data["DISK_MDOT"]=(data["L_x"]/c.c.cgs/c.c.cgs/efficiency).value   #THe disk massloss rate is only used to set the initial temperature

#Finally, lets set up the grid - confusingly, rmin and rmax are scaled by the UNIT_LENGTH, disk truncation isn't

data["R_MIN"]=0.05*R_IC/UNIT_LENGTH
data["R_MAX"]=20.*R_IC/UNIT_LENGTH
data["N_R"]=200
data["DISK_TRUNC_RAD"]=2*R_IC

data["T_MIN"]=np.radians(0.0)
data["T_MAX"]=np.radians(90.0)
data["N_T"]=100

#Now we work out the matching luminosity for the python simulation.

nu1=((13.6*u.eV).to(u.erg)/c.h.cgs).value
nu2=((2000*u.eV).to(u.erg)/c.h.cgs).value
nu3=((10000*u.eV).to(u.erg)/c.h.cgs).value
numax=(data["T_x"]*u.K*c.k_B.cgs/c.h.cgs).value*100.
L_x_test=quad(pps.brem,nu1,np.max([nu3,numax]),args=(data["T_x"],data["BREM_ALPHA"]))
const=data["L_x"]/L_x_test[0]
data["L_2_10"]=const*quad(pps.brem,nu2,nu3,args=(data["T_x"],data["BREM_ALPHA"]))[0]
print "2-10 keV luminosity=",data["L_2_10"]

data["NPHOT"]=1e7







	


t0=10000.0  #The run time for the initial zeus run - the first run is to produce a starting geometry
dt=1000.0   #
den_tol=0.5 #We ask Zeus to log cells whose density has changed by 50% or more (can be a *LOT* more)
nden=0.1    #The percentage of cells that can change before we call python again


python_ver="~/python/bin/py82k"

istart=0


if t0==0.0:
	print "We need to run for at least one second, dummy"
	t0=1.0

py_cycles=3




py_cycles=py_cycles+istart*2


out=open("pluto_py_logfile",'w',0)

out.write("Starting run"+"\n")
#out.write("zeus_ver="+zeus_ver+"\n")

for i in range(istart,10000):  #We will permit up to 500 calls to python (this is a lot)
	out.write("STARTING CYCLE "+str(i)+"\n")
	print ("STARTING CYCLE "+str(i)+"\n")
	
	pps.pluto_input_file(t0+float(i)*dt,data,rad_force)
	out.write("Running for time="+str(t0+float(i)*dt)+"\n")
	if i==0:   #This is the first step - 
		out.write("Creating first zeus_file"+"\n")
		cmdline="mpirun -n "+str(nproc_pl)+" ./pluto >"+"%08d"%i+"_pluto_log"
	else:
		out.write("generating restart zeus run \n")      #This should be the name of the restart file
		cmdline="mpirun -n "+str(nproc_pl)+" ./pluto -restart "+str(ifile)+" > "+"%08d"%i+"_pluto_log"
		
	out.write("Executing pluto with command line "+cmdline+"\n")
	subprocess.call(cmdline,shell=True)    #Call zeus
	out.write("Finished pluto run"+"\n")
	cmdline="tail -1 dbl.out"   
	out.write(cmdline+"\n")
	proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) #This mess gets the last hdffile	
	ifile=int(proc.stdout.read().split()[0])
	pps.pluto2py(ifile)   #We now make a python input file
	root="%08d"%(ifile)
	pps.python_input_file(root+".pluto",data,py_cycles)  #This generate a python parameter file
	cmdline="cp "+root+".pluto"+".pf input.pf"   #Copy the python file to a generaic name so windsave files persist
	out.write(cmdline+"\n")
	print (cmdline+"\n")

	subprocess.check_call(cmdline,shell=True)
	if py_cycles==3: #This is the first time thruogh - so no restart""
		cmdline="mpirun -n "+str(nproc_py)+" "+python_ver+" -z  input.pf > "+root+".python_log"  #We now run python
	else:
		cmdline="mpirun -n "+str(nproc_py)+" "+python_ver+" -z -r  input.pf > "+root+".python_log"  #We now run python
	out.write("Running python"+"\n") 
	print("Running python"+"\n") 	
	out.write(cmdline+"\n")
	print(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)   #Well, here is the actual call
	cmdline="cp py_heatcool.dat "+root+"_py_heatcool.dat"  
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
	py_cycles=py_cycles+2
#	now make a prefactors file
	out.write ("Makeing a prefactor file useing "+str(ifile)+" dbl file")
	pps.pre_calc(ifile)	
	cmdline="cp prefactors.dat "+root+"_prefactors.dat"  
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)
	out.write("FINISHED CYCLE"+"\n")
	
out.close()

