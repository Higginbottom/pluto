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

UNIT_TIME=UNIT_LENGTH/UNIT_VELOCITY

data={} #Set up the dictionary that will contain all of the information for the run


data["rad_force"]=1  #Including rad force? 1=yes 0=no
data["python_ver"]="84f"  #the version of python to use
data['debug']=0

data["geometry"]="spherical_1d"

data["nproc_py"]=20  #The number of cores to use for python - 256 is good!
data["nproc_pl"]=1  #The number of cores to use for pluto
data["nproc_cak"]=20  #The number of cores to use for cak


t0=50000.  #The run time for the initial pluto run - the first run is to produce a starting geometry
cycles=0
dt=-1.   #The time between calls to pluto - negative will attempt to do a fixed number of hydro timesteps
dcycles1=20 #The initial number of cycles between calls to python

ncycles1=100 #The number of cycles to use the inital cycle step for - this is to allow a new rough state to be achieved
dcycles2=1 #The number of cycles between calls to python once the initial steps are finished.


istart=0
if t0==0.0:
    print ("We need to run for at least one second, dummy")
    t0=1.0
init_py_cycles=2
py_cycles=3


#Set the parameters of the run here


data["NPHOT"]=100000

#First set the temperature of the radiation and the central source mass
data["T_star"]=42000.
data["T_iso"]=4.2e4
data["V_0"]=1e4
data["M_rad"]=1000.
data["Lum_x"]=0.0
data["Lum_UV"]=3.83e39
data["CENT_MASS"]=1.04e35
data["RHO_0"]=1.00e-11
data["R_0"]=1.317e12


#Set up the grid 


data["R_MIN"]=1.317e12
data["R_MAX"]=1.317e13
data["N_R"]=256
data["STRETCH_GRID"]=False
data["N_STRETCH"]=40
data["STRETCH_RMAX"]=300.0e10

pps.data_complete(data)


if len(sys.argv)<2:
    istart=0
    time=t0
else:
    istart=int(sys.argv[1])
    if istart>0:
        root="%08d"%(istart)
        directory="cycle"+root
        print ("We will be trying to restart using files in "+directory)
        try:
            subprocess.check_call("cp "+directory+"/* .",shell=True)
        except:
            print ("Cannot restart")
            exit(0)
        print ("Last run finished at ",pp.pload(istart).SimTime)
        if dt>0:
            time=pp.pload(istart).SimTime+dt
            cycles=-1
        else:
            time=1e99
            cmdline="tail -1 dbl.out"
            if pp.pload(istart).NStep<ncycles1:
                cycles=int(subprocess.check_output(cmdline,shell=True).split()[3])+dcycles1
            else:
                cycles=int(subprocess.check_output(cmdline,shell=True).split()[3])+dcycles2
                   
            
    else:
        istart=0
        time=t0

for i in range(istart,1):
    
    ifile=i+1
    root="%08d"%(ifile)
    print ("Making a new pluto file - time="+str(time)+" cycles="+str(cycles))
    pps.pluto_input_file(time,cycles,data)
    print ("Made pluto input file")
    
    if i==0:
        cmdline="./pluto > "+root+".pluto_log"
    else:
        if cycles==0:
            cmdline="./pluto -restart "+str(i)+" > "+root+".pluto_log"
        else:
            cmdline="./pluto -restart "+str(i)+" -maxsteps "+str(cycles+1)+" > "+root+".pluto_log"

    print ("Running pluto run")
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)
    print ("Finished pluto run")
    



    if dt>0:
        time=time+dt
        cycles=0
    else:
        cmdline="tail -1 dbl.out"
        if pp.pload(ifile).NStep<ncycles1:
            cycles=int(subprocess.check_output(cmdline,shell=True).split()[3])+dcycles1
        else:
            cycles=int(subprocess.check_output(cmdline,shell=True).split()[3])+dcycles2
        
        time=1e99


    root="%08d"%(ifile)
    dbl="data."+"%04d"%(ifile)+".dbl"
    directory="cycle"+root
    
       
    py_model_file=pps.pluto2py_1d(ifile,data)    #Make a python model file from pluto output
       
    pps.python_input_file_stellar_wind(py_model_file,data,cycles=init_py_cycles+i*py_cycles) #make a python parameter file
    subprocess.check_call("cp "+root+".pluto"+".pf input.pf",shell=True)  #Copy the python file to a generic name so windsave files persist
    print ("cycle number",i)
    
    if i>0:
        print ("restarting python")
        cmdline="modify_wind"+data["python_ver"]+" -model_file "+py_model_file+" input > "+root+".mod_wind_log" #Paint the new densities over old windsave
        print (cmdline)
        subprocess.check_call(cmdline,shell=True)   
        subprocess.check_call("cp new.wind_save input.wind_save",shell=True) 
        cmdline="mpirun -n "+str(data["nproc_py"])+" py"+data["python_ver"]+" -f -r input.pf > "+root+".python_log"        
    else:   
        cmdline="mpirun -n "+str(data["nproc_py"])+" py"+data["python_ver"]+" -f  input.pf > "+root+".python_log"
    print ("Running python")
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)
    print ("Finished python")

    cmdline="rad_hydro_files"+data["python_ver"]+" input > rad_hydro_files_output" 
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)
        
    cmdline="mpirun -n "+str(data["nproc_cak"])+" ./cak > cak_output" 
    print (cmdline)
    subprocess.check_call(cmdline,shell=True) 
    print ("Making acceleration files")
       
    pps.driving_calc(ifile,data)
    print ("Made acceleration files")

 
    if (data['debug']!=1):
        subprocess.check_call("rm -rf cycle"+"%08d"%(ifile-2),shell=True)        
    
    
    try:
        subprocess.check_call("mkdir "+directory,shell=True)
    except:
        subprocess.check_call("rm -rf "+directory+"_old",shell=True)        
        subprocess.check_call("mv "+directory+" "+directory+"_old",shell=True)
        subprocess.check_call("mkdir "+directory,shell=True)
    print ("Moving files to restart directory")

    subprocess.check_call("cp dbl.out "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp pluto.ini "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp restart.out "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp "+dbl+" "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp py_*.dat "+directory,shell=True) #Copy the rad_hydro output files to storage directory
    subprocess.check_call("cp M_data.dat "+directory,shell=True) #Copy the CAK output files to storage directory
    subprocess.check_call("cp input.wind_save "+directory,shell=True) #Copy the CAK output files to storage directory
    subprocess.check_call("cp "+py_model_file+" "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp "+root+".pluto"+".pf "+directory,shell=True)  #And also copy the python file to storage directory
    subprocess.check_call("cp "+root+".pluto_log "+directory,shell=True)  #And also copy the python file to storage directory
    subprocess.check_call("cp py_accelerations.dat "+directory,shell=True) #Copy the CAK output files to storage directory
    subprocess.check_call("cp acceleration_limits.dat "+directory,shell=True) #Copy the CAK output files to storage directory


    print ("Moved files to restart directory")



