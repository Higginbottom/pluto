#!/usr/bin/env python

import subprocess,sys
import glob
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO as pp
import numpy as np
import python_pluto_SW_sub as pps
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



try:
    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=pps.get_units()
except:
    print("Unable to open definitions.h file - big problems")
    exit()



data={} #Set up the dictionary that will contain all of the information for the run


data["rad_force"]=1  #Including rad force? 1=yes 0=no
data["python_ver"]="py84c"  #the version of python to use

data["nproc_py"]=3  #The number of cores to use for python - 256 is good!
data["nproc_pl"]=3  #The number of cores to use for pluto
data["nproc_cak"]=3  #The number of cores to use for cak


t0=1000.  #The run time for the initial pluto run - the first run is to produce a starting geometry
dt=10.0   #The time between calls to pluto
istart=0
if t0==0.0:
    print ("We need to run for at least one second, dummy")
    t0=1.0
init_py_cycles=3


#Set the parameters of the run here


data["NPHOT"]=10000

#First set the temperature of the radiation and the central source mass
data["T_star"]=50000.
data["T_iso"]=4e5
data["V_0"]=1e6
data["M_rad"]=200.
data["Lum_x"]=0.0
data["Lum_UV"]=3.77e39
data["CENT_MASS"]=60.*c.M_sun.cgs.value
data["RHO_0"]=2.00e-11
data["R_0"]=9.6e11


#Set up the grid 


data["R_MIN"]=96.0e10
data["R_MAX"]=100.0e10
data["N_R"]=100.
data["STRETCH_GRID"]=True
data["N_STRETCH"]=100
data["STRETCH_RMAX"]=10000.0e10

pps.data_complete(data)






for i in range(100):
    root="%08d"%(i)
    
    if i==0:
        pps.pluto_input_file(t0,data)
        cmdline="./pluto > "+root+".pluto_log"
        print ("Running initial pluto run")
        print (cmdline)
        subprocess.check_call(cmdline,shell=True)
        print ("Finished initial pluto run - now entering loop")
    else:
        pps.pluto_input_file(t0+float(i*dt),data)
        cmdline="./pluto -restart "+str(i)+" > "+root+".pluto_log"
        print ("Running pluto run")
        print (cmdline)
        subprocess.check_call(cmdline,shell=True)
        print ("Finished pluto run")

    ifile=i+1
    root="%08d"%(ifile)

    fname=pps.pluto2py_1d(ifile)    #Make a python model file from pluto output
    
    
    pps.python_input_file_stellar_wind(fname,data,cycles=3) #make a python parameter file
    
    cmdline="cp "+root+".pluto"+".pf input.pf"   #Copy the python file to a generaic name so windsave files persist
    subprocess.check_call(cmdline,shell=True)
    
    
    
    cmdline="mpirun -n "+str(data["nproc_py"])+" "+data["python_ver"]+" -f  input.pf > "+root+".python_log"
    print ("Running python")
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)
    print ("Finished python")

    cmdline="rad_hydro_files input > rad_hydro_files_output" 
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)

    cmdline="mpirun -n "+str(data["nproc_cak"])+" ./cak > cak_output" 
    print (cmdline)
    subprocess.check_call(cmdline,shell=True) 


    pps.driving_calc(ifile)





