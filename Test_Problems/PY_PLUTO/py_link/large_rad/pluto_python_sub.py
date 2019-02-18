#!/usr/bin/env python

import subprocess
import glob
from astropy.io import ascii
from astropy.table import Table
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad

import pyPLUTO as pp
import numpy as np


def get_units(fname='definitions.h'):
	inp=open('definitions.h','ro')
	for line in inp.readlines():
		data=line.split()
		if len(data)>1:
			if data[1]=='UNIT_DENSITY':
				UNIT_DENSITY=float(data[2])
			elif data[1]=='UNIT_LENGTH':
				UNIT_LENGTH=float(data[2])
			elif data[1]=='UNIT_VELOCITY':
				UNIT_VELOCITY=float(data[2])
	inp.close()
	return(UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY)
	
	


def brem(freq,T_x=5.6e7,alpha=0.0):
	return freq**alpha*np.exp((-1.*c.h.cgs*freq/c.k_B.cgs/T_x).value)

def pluto_input_file(tlim,data,radforce=0):
	output=open('pluto.ini','w')
	output.write("[Grid]\n")
	output.write("\n")
	output.write("X1-grid 1 "+str(data["R_MIN"])+" "+str(data["N_R"])+" r "+str(data["R_MAX"])+" 1.02\n")
	output.write("X2-grid 1 "+str(data["T_MIN"])+" "+str(data["N_T"])+" r "+str(data["T_MAX"])+" 0.95\n")
	output.write("X3-grid 1    0.0    1      u    1.0\n")
	output.write("\n")
	output.write("[Chombo Refinement]\n")
	output.write("\n")
	output.write("Levels           4\n")
	output.write("Ref_ratio        2 2 2 2 2\n") 
	output.write("Regrid_interval  2 2 2 2 \n")
	output.write("Refine_thresh    0.3\n")
	output.write("Tag_buffer_size  3\n")
	output.write("Block_factor     8\n")
	output.write("Max_grid_size    64\n")
	output.write("Fill_ratio       0.75\n")
	output.write("\n")
	output.write("[Time]\n")
	output.write("\n")
	output.write("CFL              0.4\n")
	output.write("CFL_max_var      1.1\n")
	output.write("tstop            "+str(tlim)+"\n")
	output.write("first_dt         1e-4\n")
	output.write("\n")
	output.write("[Solver]\n")
	output.write("\n")
	output.write("Solver         tvdlf\n")
	output.write("\n")
	output.write("[Boundary]\n")
	output.write("\n")
	output.write("X1-beg        outflow\n")
	output.write("X1-end        outflow\n")
	output.write("X2-beg        axisymmetric\n")
	output.write("X2-end        reflective\n")
	output.write("X3-beg        outflow\n")
	output.write("X3-end        outflow\n")
	output.write("\n")
	output.write("[Static Grid Output]\n")
	output.write("\n")
	if radforce==1:
		output.write("uservar    20    XI T ch cc lc bc xh ch_pre cc_pre lc_pre bc_pre xh_pre ne nh g1 g2 g3 g1_pre g2_pre g3_pre\n")
	else:
		output.write("uservar    14    XI T ch cc lc bc xh ch_pre cc_pre lc_pre bc_pre xh_pre ne nh\n")
	output.write("dbl        1000000000000   -1   single_file\n")
	output.write("flt       -1.0  -1   single_file\n")
	output.write("vtk       -1.0  -1   single_file\n")
	output.write("dbl.h5    -1.0  -1\n")
	output.write("flt.h5    -1.0  -1\n")
	output.write("tab       -1.0  -1   \n")
	output.write("ppm       -1.0  -1   \n")
	output.write("png       -1.0  -1\n")
	output.write("log        1000\n")
	output.write("analysis  -1.0  -1\n")
	output.write("\n")
	output.write("[Chombo HDF5 output]\n")
	output.write("\n")
	output.write("Checkpoint_interval  -1.0  0\n")
	output.write("Plot_interval         1.0  0 \n")
	output.write("\n")
	output.write("[Parameters]\n")
	output.write("\n")
	output.write("RHO_0                       "+str(data["RHO_0"])+"\n")
	output.write("RHO_ALPHA                   "+str(data["RHO_ALPHA"])+"\n")
	output.write("R_0                         "+str(data["R_0"])+"\n")
	output.write("CENT_MASS                   "+str(data["CENT_MASS"])+"\n")
	output.write("DISK_MDOT                   "+str(data["DISK_MDOT"])+"\n")
	output.write("CISO                        1e10  \n")
	output.write("L_x                         "+str(data["L_x"])+"\n")
	output.write("T_x                         "+str(data["T_x"])+"\n")
	output.write("DISK_TRUNC_RAD              "+str(data["DISK_TRUNC_RAD"])+"\n")
	output.write("MU                          "+str(data["MU"])+"\n")
	output.close()
	return
	
	
	


def python_input_file(fname,data,cycles=2):
	output=open(fname+".pf",'w')
	output.write("System_type(star,binary,agn,previous)            agn\n")
	output.write("\n")
	output.write("### Parameters for the Central Object\n")
	output.write("Central_object.mass(msol)                  7.0\n")
	output.write("Central_object.radius(cm)                  7e+08\n")
	output.write("\n")
	output.write("### Parameters for the Disk (if there is one)\n")
	output.write("\n")
	output.write("Disk.type(none,flat,vertically.extended)       none\n")
	output.write("\n")
	output.write("### Parameters for BL or AGN\n")
	output.write("\n")
	output.write("QSO_BH_radiation(yes,no)     yes\n")
	output.write("Rad_type_for_agn(0=bb,1=models,3=power_law,4=cloudy_table,5=bremsstrahlung)_to_make_wind   5\n")
	output.write("lum_agn(ergs/s) "+str(data["L_2_10"])+"\n")
	output.write("AGN.bremsstrahlung_temp(K) "+str(data["T_x"])+"\n")
	output.write("AGN.bremsstrahlung_alpha() " +str(data["BREM_ALPHA"])+"\n")
	output.write("AGN.geometry_for_pl_source(sphere,lamp_post) sphere\n")
	output.write("\n")
	output.write("### Parameters descibing the various winds or coronae in the system\n")
	output.write("\n")
	output.write("Wind_radiation(yes,no) no\n")
	output.write("Wind.number_of_components  1\n")
	output.write("Wind_type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)  hydro \n")
	output.write("Wind.coord_system(spherical,cylindrical,polar,cyl_var)  spherical\n")
	output.write("Wind.dim.in.x_or_r.direction               30\n")
	output.write("Wind.dim.in.z_or_theta.direction           30\n")
	output.write("\n")
	output.write("### Parameters associated with photon number, cycles,ionization and radiative transfer options\n")
	output.write("\n")
	output.write("Photons_per_cycle        "+str(data["NPHOT"])+"\n")
	output.write("Ionization_cycles        "+str(cycles)+"\n")
	output.write("Spectrum_cycles          0\n")
	output.write("Wind_ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
	output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
	output.write("Atomic_data  data/standard80\n")
	output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    reflect\n")
	output.write("Thermal_balance_options(0=everything.on,1=no.adiabatic)   1\n")
	output.write("\n")
	output.write("### Parameters for Domain 0\n")
	output.write("\n")
	output.write("Hydro.file "+fname+"\n")
	output.write("Hydro.thetamax(degrees:negative_means_no_maximum)  -1\n")
	output.write("Wind.t.init                                40000\n")
	output.write("Wind.filling_factor(1=smooth,<1=clumped)   1\n")
	output.write("\n")
	output.write("### Parameters for Reverberation Modeling (if needed)\n")
	output.write("\n")	
	output.write("Reverb.type(0=off,1=photon,2=wind,3=matom)   0\n")
	output.write("\n")	
	output.write("### Other parameters\n")
	output.write("\n")	
	output.write("Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)  logarithmic\n")
	output.write("Photon_sampling.nbands                     10\n")
	output.write("Photon_sampling.low_energy_limit(eV)       1.03333\n")
	output.write("Photon_sampling.high_energy_limit(eV)      50000\n")
	
	output.close()
	return
		
	
def pluto2py(ifile):

	D=pp.pload(ifile)

	# We need the definitions file - so we know the conversion factors.

	UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')

	# Open an output file 

	fname="%08d"%ifile+".pluto"

	# Preamble

	out=open(fname,'w')
	out.write("# This is a file generated by hydro_to_python\n")
	out.write("# We can put any number of comments in behind # signs\n")
	out.write("# By default, the order of coordinates are \n")
	out.write("#				r, theta phi for spherical polars\n")
	out.write("# 		                x,y,z        for carteisan\n")
	out.write("#                         or w, z, phi    for cylindrical\n")


	titles=[]
	titles=titles+["ir","r_cent","r_edge"]
	titles=titles+["itheta","theta_cent","theta_edge"]
	titles=titles+["v_r","v_theta","v_phi","density","temperature"]

	r_edge=[]
	r_ratio=(D.x1[2]-D.x1[1])/(D.x1[1]-D.x1[0])
	dr=(D.x1[1]-D.x1[0])/(0.5*(1.0+r_ratio))
	r_edge.append(D.x1[0]-0.5*dr)
	for i in range(len(D.x1)-1):
		r_edge.append(r_edge[-1]+dr)
		dr=dr*r_ratio
	
	
	r_edge=np.array(r_edge)	

	theta_edge=[]
	theta_ratio=(D.x2[2]-D.x2[1])/(D.x2[1]-D.x2[0])
	dtheta=(D.x2[1]-D.x2[0])/(0.5*(1.0+theta_ratio))
	theta_min=D.x2[0]-0.5*dtheta
	if theta_min<0.0:
		theta_min=0.0
	theta_edge.append(theta_min)
	for i in range(len(D.x2)-1):
		theta_edge.append(theta_edge[-1]+dtheta)
		dtheta=dtheta*theta_ratio
	if (theta_edge[-1]+(D.x2[-1]-theta_edge[-1])*2.0)>(np.pi/2.0):
		D.x2[-1]=(theta_edge[-1]+(np.pi/2.0))/2.0

	theta_edge=np.array(theta_edge)	

	col0=np.array([])
	col1=np.array([])
	col2=np.array([])
	col3=np.array([])
	col4=np.array([])
	col5=np.array([])
	col6=np.array([])
	col7=np.array([])
	col8=np.array([])
	col9=np.array([])
	col10=np.array([])

	fmt='%013.6e'

	#This next line defines formats for the output variables. This is set in a dictionary
	fmts={	'ir':'%03i',	
		'r_cent':fmt,
		'r_edge':fmt,
		'itheta':'%i',	
		'theta_cent':fmt,
		'theta_edge':fmt,
		'v_r':fmt,
		'v_theta':fmt,
		'v_phi':fmt,
		'density':fmt,
		'temperature':fmt}

	for j in range(len(D.x2)):
		col0=np.append(col0,np.arange(len(D.x1)))
		col1=np.append(col1,D.x1*UNIT_LENGTH)
		col2=np.append(col2,r_edge*UNIT_LENGTH)
		col3=np.append(col3,np.ones(len(D.x1))*j)
		col4=np.append(col4,np.ones(len(D.x1))*D.x2[j])
		col5=np.append(col5,np.ones(len(D.x1))*theta_edge[j])
		col6=np.append(col6,np.transpose(D.vx1)[j]*UNIT_VELOCITY)
		col7=np.append(col7,np.transpose(D.vx2)[j]*UNIT_VELOCITY)
		col8=np.append(col8,np.transpose(D.vx3)[j]*UNIT_VELOCITY)
		col9=np.append(col9,np.transpose(D.rho)[j]*UNIT_DENSITY)
		col10=np.append(col10,np.transpose(D.T)[j])

	out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10],names=titles)
	ascii.write(out_dat,out,formats=fmts)
	out.close()
	return


def pre_calc(ifile,radforce=0):
	max_change=0.9
	max_accel_change=0.9

	heatcool=ascii.read("%08d"%(ifile)+"_py_heatcool.dat")
	D=pp.pload(ifile)

	# We need the definitions file - so we know the conversion factors.

	UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')
	UNIT_ACCELERATION=UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH	

	comp_h_pre=[]
	comp_c_pre=[]
	xray_h_pre=[]
	brem_c_pre=[]
	line_c_pre=[]
	
	g1_pre=[]
	g2_pre=[]
	g3_pre=[]

	odd=0.0
	
	itest=19900

	for i in range(len(heatcool["rho"])):
		if (heatcool["rho"][i]/(D.rho[heatcool["i"][i]][heatcool["j"][i]]*UNIT_DENSITY))-1.>1e-6:
			odd=odd+1
		nenh=D.ne[heatcool["i"][i]][heatcool["j"][i]]*D.nh[heatcool["i"][i]][heatcool["j"][i]]
		nhnh=D.nh[heatcool["i"][i]][heatcool["j"][i]]*D.nh[heatcool["i"][i]][heatcool["j"][i]]
		
		ideal_prefactor=(heatcool["heat_comp"][i]/(D.ch[heatcool["i"][i]][heatcool["j"][i]]*nenh))
		change=ideal_prefactor/D.ch_pre[heatcool["i"][i]][heatcool["j"][i]]
		if change<max_change:
			change=max_change
		elif change>(1./max_change):
			change=(1./max_change)
		comp_h_pre.append(change*D.ch_pre[heatcool["i"][i]][heatcool["j"][i]])
			
		ideal_prefactor=(heatcool["cool_comp"][i]/(D.cc[heatcool["i"][i]][heatcool["j"][i]]*nenh))
		change=ideal_prefactor/D.cc_pre[heatcool["i"][i]][heatcool["j"][i]]
		if change<max_change:
			change=max_change
		elif change>(1./max_change):
			change=(1./max_change)
		comp_c_pre.append(change*D.cc_pre[heatcool["i"][i]][heatcool["j"][i]])
	
		ideal_prefactor=(heatcool["cool_lines"][i]/(D.lc[heatcool["i"][i]][heatcool["j"][i]]*nenh))
		change=ideal_prefactor/D.lc_pre[heatcool["i"][i]][heatcool["j"][i]]
		if change<max_change:
			change=max_change
		elif change>(1./max_change):
			change=(1./max_change)
		line_c_pre.append(change*D.lc_pre[heatcool["i"][i]][heatcool["j"][i]])
		
		ideal_prefactor=(heatcool["cool_ff"][i]/(D.bc[heatcool["i"][i]][heatcool["j"][i]]*nenh))
		change=ideal_prefactor/D.bc_pre[heatcool["i"][i]][heatcool["j"][i]]
		if change<max_change:
			change=max_change
		elif change>(1./max_change):
			change=(1./max_change)
		brem_c_pre.append(change*D.bc_pre[heatcool["i"][i]][heatcool["j"][i]])
	
		ideal_prefactor=(heatcool["heat_xray"][i]/(D.xh[heatcool["i"][i]][heatcool["j"][i]]*nenh))
		change=ideal_prefactor/D.xh_pre[heatcool["i"][i]][heatcool["j"][i]]
		if change<max_change:
			change=max_change
		elif change>(1./max_change):
			change=(1./max_change)
		xray_h_pre.append(change*D.xh_pre[heatcool["i"][i]][heatcool["j"][i]])
		
		
		if radforce:
			ideal_prefactor=(heatcool["rad_f_w"][i]/heatcool["rho"][i]/heatcool["vol"][i])/(D.g1[heatcool["i"][i]][heatcool["j"][i]]*UNIT_ACCELERATION)
			change=ideal_prefactor/D.g1_pre[heatcool["i"][i]][heatcool["j"][i]]
			if change<max_accel_change:
				change=max_accel_change
			elif change>(1./max_accel_change):
				change=(1./max_accel_change)
			g1_pre.append(change*D.g1_pre[heatcool["i"][i]][heatcool["j"][i]])
			
			g2_pre.append(1.0)
			
			ideal_prefactor=(heatcool["rad_f_z"][i]/heatcool["rho"][i]/heatcool["vol"][i])/(D.g3[heatcool["i"][i]][heatcool["j"][i]]*UNIT_ACCELERATION)
			change=ideal_prefactor/D.g3_pre[heatcool["i"][i]][heatcool["j"][i]]
			if change<max_accel_change:
				change=max_accel_change
			elif change>(1./max_accel_change):
				change=(1./max_accel_change)
			g3_pre.append(change*D.g3_pre[heatcool["i"][i]][heatcool["j"][i]])	
		else:		
			g1_pre.append(1.0)
			g2_pre.append(1.0)
			g3_pre.append(1.0)
			
	
	fmt='%013.6e'

	#This next line defines formats for the output variables. This is set in a dictionary
	fmts={	'ir':'%03i',
		'rcent':fmt,
		'itheta':'%03i',
		'thetacent':fmt,	
		'rho':fmt,
		'comp_h_pre':fmt,
		'comp_c_pre':fmt,
		'xray_h_pre':fmt,
		'line_c_pre':fmt,
		'brem_c_pre':fmt,
		'g1_pre':fmt,
		'g2_pre':fmt,
		'g3_pre':fmt,
		}	

	titles=[]
	titles=titles+["ir","rcent","itheta","thetacent","rho"]
	titles=titles+["comp_h_pre","comp_c_pre","xray_h_pre","brem_c_pre","line_c_pre"]
	titles=titles+["g1_pre","g2_pre","g3_pre"]	
	
	col0=heatcool["i"]
	col1=heatcool["rcen"]
	col2=heatcool["j"]
	col3=heatcool["thetacen"]
	col4=heatcool["rho"]
	col5=comp_h_pre
	col6=comp_c_pre
	col7=xray_h_pre
	col8=brem_c_pre
	col9=line_c_pre
	col10=g1_pre
	col11=g2_pre
	col12=g3_pre

	out=open("prefactors.dat",'w')

	out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12],names=titles)
	ascii.write(out_dat,out,formats=fmts)
	return(odd)