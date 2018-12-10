#!/usr/bin/env python 

import subprocess
import glob



def python_input_file(fname,cycles=2):
	output=open(fname+".pf",'w')
	output.write("System_type(0=star,1=binary,2=agn)                    2\n")
	output.write("Wind_type                 3\n")	
	output.write("Coord.system(0=spherical,1=cylindrical,2=spherical_polar,3=cyl_var)                    1\n")
	output.write("Wind.dim.in.x_or_r.direction                     30\n")
	output.write("Wind.dim.in.z_or_theta.direction                   30\n")
	output.write("Number.of.wind.components 1\n") 
	output.write("disk.type(0=no.disk,1=standard.flat.disk,2=vertically.extended.disk) 0\n") 
	output.write("Atomic_data                         data/standard79_coll\n")
	output.write("write_atomicdata(0=no,1=yes)               0\n")        
	output.write("photons_per_cycle                            10000000\n")
	output.write("Ionization_cycles                                "+str(cycles)+"\n")
	output.write("spectrum_cycles                                   0\n")
	output.write("adjust_grid(0=no,1=yes)								0\n")
	output.write("Wind_ionization 9\n")
	output.write("Line_transfer 3\n")
	output.write("Thermal_balance_options(0=everything.on,1=no.adiabatic)                    1\n")
	output.write("Disk_radiation(y=1)                               0\n")
	output.write("Wind_radiation(y=1)                               1\n")
	output.write("QSO_BH_radiation                               1\n")
	output.write("Rad_type_for_disk(0=bb,1=models)_to_make_wind     0\n")
	output.write("Rad_type_for_agn(0=bb,1=models,3=power_law,4=cloudy_table)_to_make_wind)  5\n")
	output.write("mstar(msol)                                     7\n")
	output.write("rstar(cm)                                     7e+08\n")
	output.write("tstar                                         40000\n")
	output.write("lum_agn(ergs/s) 1.7698e37\n")
	output.write("agn_bremsstrahlung_temp(K) 5.6e7\n")
	output.write("agn_bremsstrahlung_alpha 0.0\n")
	output.write("geometry_for_pl_source 0\n")
	output.write("agn_power_law_index 							0.0\n")
	output.write("agn_power_law_cutoff (0)						0\n")
	output.write("Torus(0=no,1=yes)								0\n")
	output.write("disk.mdot(msol/yr)                            1e-08\n")
	output.write("Disk.illumination.treatment 0\n")
	output.write("Disk.temperature.profile(0=standard;1=readin)                    0\n")
	output.write("disk.radmax(cm)                             2.4e+10\n")
	output.write("wind.radmax(cm)                               1e+11\n")
	output.write("wind.t.init                                   40000\n")
	output.write("hydro_file "+fname+".zeu\n")
	output.write("Hydro_thetamax(degrees)                        -1\n")
	output.write("filling_factor(1=smooth,<1=clumped)                    1\n")
	output.write("Rad_type_for_agn(3=power_law,4=cloudy_table)_in_final_spectrum 3\n")
	output.write("Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum                    0\n")
	output.write("spectrum_wavemin                               1450\n")
	output.write("spectrum_wavemax                               1650\n")
	output.write("no_observers                                      4\n")
	output.write("angle(0=pole)                                    10\n")
	output.write("angle(0=pole)                                    30\n")
	output.write("angle(0=pole)                                    60\n")
	output.write("angle(0=pole)                                    80\n")
	output.write("live.or.die(0).or.extract(anything_else)                    1\n")
	output.write("spec.type(flambda(1),fnu(2),basic(other)                    1\n")
	output.write("Use.standard.care.factors(1=yes)						1\n")
	output.write("reverb.type 0\n")
	output.write("Photon.sampling.approach           8\n")
	output.write("Num.of.frequency.bands(5) 10\n")
	output.write("Lowest_energy_to_be_considered(eV) 1.03333\n")
	output.write("Highest_energy_to_be_considered(eV) 50000 \n")
	output.write("Extra.diagnostics(0=no,1=yes)   1\n")
	output.write("keep_ioncycle_windsaves()   1\n")
	output.close()
	return


def zeus_input_file(resfile,restart,tlim,den_tol=0.5,nden=0.1,hc_lim=0.9,runname='aa'):
	output=open("z2dinput","w")
	output.write("&rescon  irestart="+str(restart)+",resfile='"+resfile+"',dtdump=72000.0,id='"+runname+"' /\n")
	output.write("&pcon nlim=10000000, tlim="+str(tlim)+" /\n")
	output.write("&hycon qcon=2.0,courno=0.5 /\n")
	output.write("&ggen1 nbl=80,x1min=2.4123382003e10,x1max=96.4935280142e10 ,igrid=1,x1rat=1.02,lgrid=.true.  /\n")
	output.write("&ggen2 nbl=100,x2min=0.0,x2max=1.5707964,igrid=1,x2rat=0.95,lgrid=.true.  /\n")
        output.write("&iib niib=205*1, miib=205*1  /\n")
        output.write("&oib noib=205*2, moib=205*2  /\n")
        output.write("&ijb nijb=205*-1,mijb=205*-1  /\n")
        output.write("&ojb nojb=205*1,mojb=205*1  /\n")
 	output.write("&grvcon ptmass=13.9e33  /\n")
 	output.write("&radcon /\n")
 	output.write("&eqos gamma=1.66666 /\n")
 	output.write("&pgen d0=1.6e-12,v10=0.,v20=0.,v30=0.,dindex=2.0,mdot=4.4e17,tx=5.6e7, lx=3.3e37, mu=0.6/\n")
 	output.write("&hc_pre comp_pre0=1.0,line_pre0=1.0,brem_pre0=1.0,xray_pre0=1.0 /\n")
	output.write("&py_hydro den_tol="+str(den_tol)+",nden="+str(nden)+",hc_lim="+str(hc_lim)+" /\n")
 	output.write("&gcon /\n")
 	output.write("&iocon dthist=1.0,dtusr=-1,dthdf=1000.0 /\n")
	output.close()
	return
	
#subprocess.check_call("touch usr",shell=True)
#subprocess.check_call("rm usr*",shell=True)
#subprocess.check_call("touch res",shell=True)
#subprocess.check_call("rm res*",shell=True)
#subprocess.check_call("touch hdf",shell=True)
#subprocess.check_call("rm hdf*",shell=True)
#subprocess.check_call("touch hst",shell=True)
#subprocess.check_call("rm hst*",shell=True)

#subprocess.check_call("touch py_heatcool.dat",shell=True)
#subprocess.check_call("rm py_heatcool.dat",shell=True)


t0=10000.0  #The run time for the initial zeus run - the first run is to produce a starting geometry
dt=1000.0
den_tol=0.5 #We ask Zeus to log cells whose density has changed by 50% or more (can be a *LOT* more)
nden=0.1    #The percentage of cells that can change before we call python again

zeus_ver="/home/nsh2m14/zeus/exe/zeus2d_lu_trunc.exe"
python_ver="/home/nsh2m14/python/bin/py81c_dev"

istart=0
tlim=1000000  #

runname='nsh_dw44'

if t0==0.0:
	print "We need to run for at least one second, dummy"
	t0=1.0

py_cycles=3


py_cycles=py_cycles+istart*2


out=open("python_logfile",'w',0)

out.write("Starting run"+"\n")
out.write("zeus_ver="+zeus_ver+"\n")

for i in range(istart,10000):  #We will permit up to 500 calls to python (this is a lot)
	out.write("STARTING CYCLE "+str(i)+"\n")
	if i==0:   #This is the first step - 
		out.write("Creating first zeus_file"+"\n")
		zeus_input_file ("",0,t0,0.5,1.1,0.9,runname)  #There is no restart file, 0 says this isnt a restart, we need to at least run for 1 second
	else:
		fname='res'+str(i-1).zfill(6)+runname  #This is a step which we expect a python run to have already taken place
		out.write("generating restart zeus run, restart file="+fname+"\n")      #This should be the name of the restart file
		zeus_input_file (fname,1,t0+dt*float(i),den_tol,nden,0.9,runname)  #We generate a new zeus input file - this is a restart
	cmdline=zeus_ver+" > zeus_"+str(i)+".out"	
	out.write("Executing zeus with command line "+cmdline+"\n")
	subprocess.call(cmdline,shell=True)    #Call zeus
	cmdline="cp z2dinput z2dinput_"+str(i)   #We will retain all the z2d input files for later
	out.write("Finished zeus run"+"\n")
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)
	cmdline="ls -tr ./hdf*"+runname+" > temp"  #Its kind of tricky to know which the last hdf file is - this approach is risky!!
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)  #but it should make a file with the most recent hdf*aa file as the last entry
	cmdline="tail -1 temp"   
	out.write(cmdline+"\n")
	proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) #This mess gets the last hdffile
	hdf_file=proc.stdout.read()[:-1]  #Thus should be the hdf file we are interested in
	cmdline="hydro_2_python.py > output "+hdf_file[2:]   #We now make a python input file
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)
	python_input_file(hdf_file[2:],py_cycles)  #This generate a python parameter file
	cmdline="cp "+hdf_file[2:]+".pf input.pf"   #Copy the python file to a generaic name so windsave files persist
	out.write(cmdline+"\n")
	
	subprocess.check_call(cmdline,shell=True)
	if py_cycles==3: #This is the first time thruogh - so no restart""
		cmdline="mpirun -n 256 "+python_ver+" -z  input.pf > "+hdf_file[2:]+".out"  #We now run python
	else:
		cmdline="mpirun -n 256 "+python_ver+" -z -r  input.pf > "+hdf_file[2:]+".out"  #We now run python
			
	out.write("Running python"+"\n") 
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)   #Well, here is the actual call
	cmdline="cp py_heatcool.dat "+hdf_file[2:]+"_py_heatcool.dat"  
	out.write(cmdline+"\n")
	subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
#	cmdline="rm "+hdf_file[2:]+".wind_save"  
#	subprocess.check_call(cmdline,shell=True)   #Remove windsave files.
	py_cycles=py_cycles+2
	out.write("FINISHED CYCLE"+"\n")
	
out.close()

