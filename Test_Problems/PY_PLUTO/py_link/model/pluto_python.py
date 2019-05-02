#!/usr/bin/env python

import subprocess,sys
import glob
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO as pp
import numpy as np
import pluto_python_sub as pps
import re

help= '''
This is pluto_python - it needs to be called with at least one option

The options are

help         - show this message
start_test   - test that the setup will run 
start        - start a new run
restart      - restart from an existing run
restart_test - test the point at which a last run finished, and suggest restart parameters
clean        - remove all runfiles - take care!!!
'''




data={} #Set up the dictionary that will contain all of the information for the run


data["rad_force"]=0  #Including rad force? 1=yes 0=no
data["python_ver"]="/Users/nsh2m14/python/bin/py83a"  #the version of python to use
data["nproc_py"]=1  #The number of cores to use for python - 256 is good!
data["nproc_pl"]=1  #The number of cores to use for pluto

t0=100.0  #The run time for the initial pluto run - the first run is to produce a starting geometry
dt=100.0   #The time between calls to pluto
istart=0
if t0==0.0:
	print "We need to run for at least one second, dummy"
	t0=1.0
init_py_cycles=3


#Set the parameters of the run here

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

data["L_x"]=3.3e37       #The luminosity from 13.6eV to infinity


#data["spectype"]='brem'
#data["BREM_ALPHA"]=0.0   #The power law for the bremstrahlung spectrum - should stay at zero unless there is a very goo reason


data["spectype"]='models'
data["model"]='model.ls'   #The power law for the bremstrahlung spectrum - should stay at zero unless there is a very goo reason


#Finally, lets set up the grid 

data["R_MIN"]=0.5*R_IC
data["R_MAX"]=20.*R_IC
data["N_R"]=20
data["DISK_TRUNC_RAD"]=2*R_IC

data["T_MIN"]=np.radians(0.0)
data["T_MAX"]=np.radians(90.0)
data["N_T"]=100

data["NPHOT"]=1e4

pps.data_complete(data) #Compute various dependant variables



if len(sys.argv)<2:
	print "We need one option to proceed"
	print help
	exit()
	
	
if sys.argv[1]=="start_test":
	print "Testing setup and generating script file"
	cmdline=python_ver+" -i input > output"
	try:
		subprocess.call(cmdline,shell=True)
		print "python runs"
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			print python_ver,"Doesn't exist"
		else:
        # Something else went wrong while trying to run `wget`
			print "Some error occured while trying to run python"
	exit()
		

elif sys.argv[1]=="restart_test":
	dbl_file_1,dbl_file_2,py_last_completed,py_last_requested=pps.get_status()
	
	if dbl_file_1!=dbl_file_2:
		print ("There is a mismatch between the dbl file in the directory and in the dbl.out file - needs checking")
	
	py_test=3+(dbl_file_1-1)*2
	print "Last requested python cycle :",py_last_requested
	print "Last completed python cycle :",py_last_completed
	print "Starting from dbl",dbl_file_1,"would normally call python with",py_test,"cycles"

	if py_last_requested==py_last_completed:
		print "I think that the sim was interrupted during a pluto run"
		py_cycles=py_last_requested+2 
		istart=np.min([dbl_file_1,dbl_file_2])		
		print "Suggest restarting with command pluto_python restart ",istart,py_cycles	
		print 3+(istart-1)*2,py_cycles			
	elif (py_last_requested-py_last_completed)==2:
		print "I think that the sim was interrupted before a python run had done one cycle"
		print "So the last pluto .dbl file is different from the densities in the last prefactor"
		py_cycles=py_last_requested	
		istart=np.min([dbl_file_1,dbl_file_2])
		print "Suggest restarting with command pluto_python restart ",istart,py_cycles		
	elif (py_last_requested-py_last_completed)==1:
		print "I think that the sim was interrupted midway through a python run"
		print "We need to restart the python run"
		py_cycles=py_last_requested	
		istart=np.min([dbl_file_1,dbl_file_2])
		print "Suggest restarting with command pluto_python restart ",istart,py_cycles
	exit()
	
	
elif sys.argv[1]=="clean":
	request= raw_input("Are you sure (Y/N)?")
	if request=="Y":
		try:
			subprocess.check_call("rm dbl*",shell=True)
			subprocess.check_call("rm data*.dbl",shell=True)
			subprocess.check_call("rm *_prefactors.dat",shell=True)		
			subprocess.check_call("rm *_py_heatcool.dat",shell=True)		
			subprocess.check_call("rm *_pluto_log",shell=True)		
			subprocess.check_call("rm *.pluto",shell=True)		
			subprocess.check_call("rm *.pluto.pf",shell=True)		
			subprocess.check_call("rm *.python_log",shell=True)		
			subprocess.check_call("rm *.wind_save",shell=True)
			subprocess.check_call("rm -r diag*",shell=True)			
			exit()
		except:
			print "Clean!"
			exit()
	else:
		exit()
		
elif sys.argv[1]=="start":
	istart=0
	if t0==0.0:
		print "We need to run for at least one second, dummy"
		t0=1.0

	py_cycles=3 #We start off with three python cycles
	istart=0 #This is a start		
		
		
elif sys.argv[1]=="restart":
	if len(sys.argv)<4:
		print "Need to include 1: dbl file to restart from 2: number of python cycles to do"
	else:
		istart=int(sys.argv[2])	
		py_cycles=int(sys.argv[3])	
		
	dbl_file_1,dbl_file_2,py_last_completed,py_last_requested=pps.get_status()
	
	py_test=3+(istart-1)*2	
	
	print py_test,py_last_requested
	
	if py_test==py_last_requested:
		print "We will be starting off by completing the last python run"
		


elif sys.argv[1]=="help" or sys.argv[1]=='h':
	print help
	exit ()
	
else:
	print "Don't understand the input - exiting"
	exit()

	








pps.loop(t0,dt,istart,py_cycles,data)
