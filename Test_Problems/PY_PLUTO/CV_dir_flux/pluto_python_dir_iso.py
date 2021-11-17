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
import time as timer

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


disk_rad_eff=0.06 #disk reprocessing efficiency
eddington_ratio=0.5


data={} #Set up the dictionary that will contain all of the information for the run

data["nproc_py"]=5
data["nproc_cak"]=5



#First, set up the grid 

data["R_MIN"]=8.7e8
data["R_MAX"]=87.e8
data["N_R"]=128

data["DISK_TRUNC_RAD"]=87.e8

data["T_MIN"]=np.radians(0.0)
data["T_MAX"]=np.radians(90.0)
data["N_T"]=96



#Set the parameters of the run here

#First set the temperature of the radiation and the central source mass

data["rad_force"]=1
data["RHO_0"]=1e-8
data["RHO_ALPHA"]=0.0
data["R_0"]=8.31e17 



data["T_x"]=160000.
data["T_star"]=40000.
data["CENT_MASS"]=0.6*c.M_sun.cgs.value
data["MU"]=0.6


L_edd=1.26e38*data["CENT_MASS"]/c.M_sun.cgs.value


data["CENT_RADIUS"]=8.7e8

data["system_type"]='star'
data["boundary_layer"]="yes"

data["L_star"]=9.05e34 #Used for the inital pluto run
data["L_BL"]=9.05e+34
data["T_BL"]=40000.

data["f_uv"]=0.9
data["f_x"]=0.1 #This is the proportion of disk flux we assume comes from the middle as X-rays

#Stuff for the python runs

data["NPHOT"]=1e4
data["PY_DISK_MDOT"]=3.14e-8
data["DISK_MDOT"]=data["PY_DISK_MDOT"]*c.M_sun.cgs.value/60/60/24/365.25
data["disk_radiation"]="yes"
data["cent_spectype"]="none" #Turns off the python central source
data["python_ver"]="85i_rh"  #the version of python to use

#Rescale the grid
    
data["R_MIN"]=data["R_MIN"]/UNIT_LENGTH
data["R_MAX"]=data["R_MAX"]/UNIT_LENGTH




t0=1.0 #The run time for the initial pluto run - the first run is to produce a starting geometry
dt=1.0   #The time between calls to pluto (in seconds)
istart=0
if t0==0.0:
    print ("We need to run for at least one second, dummy")
    t0=1.0
init_py_cycles=1
py_cycles=1




logfile=open("run_logfile","w")


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
        time=pp.pload(istart).SimTime+dt
    else:
        istart=0
        time=t0



data["k"]=0.59
data["alpha"]=-0.6



for i in range(istart,10000):
    if i>0:
        data["k"]=999
        data["alpha"]=999
        
    timer0=timer.perf_counter()
        
    root="%08d"%(i)
    logfile.write("Making a pluto input file for cycle "+str(i)+"\n")
    print("Making a pluto input file for cycle "+str(i)+"\n")
    
    pps.pluto_input_file(time,data)
    logfile.write("starting "+str(timer.perf_counter()-timer0)+"\n")
    print("starting "+str(timer.perf_counter()-timer0)+"\n")
    
    if i==0:
        cmdline="./pluto > pluto_log"
    else:
        cmdline="./pluto -restart "+str(i)+" > pluto_log"

    logfile.write("Running pluto run"+"\n")
    print("Running pluto run"+"\n")
    
    logfile.write("command line: "+cmdline+"\n")
    print("command line: "+cmdline+"\n")
    
    subprocess.check_call(cmdline,shell=True)
    logfile.write("Finished pluto run"+"\n")
    print("Finished pluto run"+"\n")
    
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
    print(str(timer.perf_counter()-timer0)+"\n")

   

    ifile=i+1   #The current pluto save file will by one further on....
    time=time+dt #The next pluto run will need to run a little longer
    root="%08d"%(ifile) #This is the root we will use for all file associated with this pluto save file
    dbl="data."+"%04d"%(ifile)+".dbl"  #The name of the dbl file 
    directory="cycle"+root #The name of the directory we will save all the files into

    logfile.write("Turning dbl file "+dbl+" into a python model file"+"\n")
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
    
    print("Turning dbl file "+dbl+" into a python model file"+"\n")
    print(str(timer.perf_counter()-timer0)+"\n")
    
    
    
    py_model_file=pps.pluto2py_rtheta(ifile)    #Make a python model file from pluto output
    logfile.write("Made a python model file called "+py_model_file+"\n")
    logfile.write(str(timer.perf_counter()-timer0)+"\n")

    print("Made a python model file called "+py_model_file+"\n")
    print(str(timer.perf_counter()-timer0)+"\n")

    logfile.write("Making a python input file"+"\n")
    print("Making a python input file"+"\n")
    
    pps.python_input_file(root,data,cycles=init_py_cycles+i*py_cycles)
    logfile.write("Successfully made a python input file"+"\n")
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
    print("Successfully made a python input file"+"\n")
    print(str(timer.perf_counter()-timer0)+"\n")

    cmdline="cp "+root+".pf input.pf"
    logfile.write("command line: "+cmdline+"\n")
    print("command line: "+cmdline+"\n")
    
    subprocess.check_call(cmdline,shell=True)  #Copy the python file to a generic name so windsave files persist  



    if i>0:
        logfile.write("restarting python"+"\n")
        logfile.write("running modify wind on the old windsave"+"\n")        
        cmdline="modify_wind"+data["python_ver"]+" -model_file "+py_model_file+" input > mod_wind_log" #Paint the new densities over old windsave
        logfile.write("command line: "+cmdline+"\n")
        subprocess.check_call(cmdline,shell=True) 
        cmdline="cp new.wind_save input.wind_save"
        logfile.write("command line: "+cmdline+"\n")          
        subprocess.check_call(cmdline,shell=True) 
        cmdline="mpirun -n "+str(data["nproc_py"])+" py"+data["python_ver"]+" -f -r -classic input.pf > python_log"        
    else:   
        cmdline="mpirun -n "+str(data["nproc_py"])+" py"+data["python_ver"]+" -f -classic input.pf > python_log"
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
        
    logfile.write("Running python"+"\n")
    logfile.write("command line: "+cmdline+"\n")
    print("Running python"+"\n")
    print("command line: "+cmdline+"\n")          
    subprocess.check_call(cmdline,shell=True)
    logfile.write("Finished python"+"\n")
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
    print("Finished python"+"\n")
    print(str(timer.perf_counter()-timer0)+"\n")

    cmdline="rad_hydro_files"+data["python_ver"]+" input > rad_hydro_files_output" 
    logfile.write("command line: "+cmdline+"\n") 
    logfile.write(str(timer.perf_counter()-timer0)+"\n")
    print("command line: "+cmdline+"\n") 
    print(str(timer.perf_counter()-timer0)+"\n")         
    subprocess.check_call(cmdline,shell=True)

    cmdline="cp py_heatcool.dat "+root+"_py_heatcool.dat"  
    subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later 

    cmdline="cp py_driving.dat "+root+"_py_driving.dat"  
    subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later 

    cmdline="cp py_pcon_data.dat "+root+"_py_pcon_data.dat"  
    subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later 


    subprocess.check_call("cp directional_flux_x_opt_thin.dat directional_flux_x.dat",shell=True)   
    subprocess.check_call("cp directional_flux_y_opt_thin.dat directional_flux_y.dat",shell=True)   
    subprocess.check_call("cp directional_flux_z_opt_thin.dat directional_flux_z.dat",shell=True)   



    cmdline="mpirun -n "+str(data["nproc_cak"])+" ./cak_v2 > cak_output"
    print ("Running CAK")
    print (cmdline)
    subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later 
    print ("Finished CAK")
    
     
    subprocess.check_call("rm -rf cycle"+"%08d"%(ifile-2),shell=True)        
    
    
    try:
        subprocess.check_call("mkdir "+directory,shell=True)
    except:
        subprocess.check_call("rm -rf "+directory+"_old",shell=True)        
        subprocess.check_call("mv "+directory+" "+directory+"_old",shell=True)
        subprocess.check_call("mkdir "+directory,shell=True)
    
    subprocess.check_call("cp dbl.out "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp pluto.ini "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp restart.out "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp "+dbl+" "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("cp py_*.dat "+directory,shell=True) #Copy the rad_hydro output files to storage directory
    subprocess.check_call("cp M_UV_data.dat "+directory,shell=True) #Copy the CAK output files to storage directory
    logfile.write("B4"+str(timer.perf_counter()-timer0)+"\n")
    
    subprocess.check_call("cp input.wind_save "+directory,shell=True) #Copy the CAK output files to storage directory
    logfile.write("AF"+str(timer.perf_counter()-timer0)+"\n")
    
    subprocess.check_call("mv "+py_model_file+" "+directory,shell=True) #Copy the model file to the storage directory
    subprocess.check_call("mv "+root+".pf "+directory,shell=True)  #And also copy the python file to storage directory
    subprocess.check_call("mv "+root+"_py_heatcool.dat "+directory,shell=True)  
    subprocess.check_call("mv "+root+"_py_driving.dat "+directory,shell=True)  
    subprocess.check_call("mv "+root+"_py_pcon_data.dat "+directory,shell=True)  
    
    #   subprocess.check_call("cp prefactors.dat "+directory,shell=True) #Copy the CAK output files to storage directory
    logfile.write("Finished tidying up"+"\n")
    logfile.write("finshed"+str(timer.perf_counter()-timer0)+"\n")
    
    
print ("Fin")