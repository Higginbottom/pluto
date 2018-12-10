#!/usr/bin/env python 

import subprocess,sys
import glob
from astropy.io import ascii
from astropy.table import Table
import pyPLUTO as pp
import numpy as np


ifile=int(sys.argv[1])

max_change=0.9
heatcool=ascii.read("py_heatcool.dat")
D=pp.pload(ifile)

# We need the definitions file - so we know the conversion factors.

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

comp_h_pre=[]
comp_c_pre=[]
xray_h_pre=[]
brem_c_pre=[]
line_c_pre=[]

odd=0.0

for i in range(len(heatcool["rho"])):
	if (heatcool["rho"][i]/(D.rho[heatcool["i"][i]][heatcool["j"][i]]*UNIT_DENSITY))-1.>1e-6:
		odd=odd+1
	nenh=D.ne[heatcool["i"][i]][heatcool["j"][i]]*D.nh[heatcool["i"][i]][heatcool["j"][i]]
	test=(heatcool["heat_comp"][i]/(D.comp_h_pre[heatcool["i"][i]][heatcool["j"][i]]*D.comp_h[heatcool["i"][i]][heatcool["j"][i]]*nenh))
	if test<max_change*D.comp_h_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=max_change*D.comp_h_pre[heatcool["i"][i]][heatcool["j"][i]]
	elif test>(1./max_change)*D.comp_h_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=(1./max_change)*D.comp_h_pre[heatcool["i"][i]][heatcool["j"][i]]
	comp_h_pre.append(test)
		
	test=(heatcool["cool_comp"][i]/(D.comp_c_pre[heatcool["i"][i]][heatcool["j"][i]]*D.comp_c[heatcool["i"][i]][heatcool["j"][i]]*nenh))
	if test<max_change*D.comp_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=max_change*D.comp_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	elif test>(1./max_change)*D.comp_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=(1./max_change)*D.comp_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	comp_c_pre.append(test)	

	test=(heatcool["cool_lines"][i]/(D.line_c_pre[heatcool["i"][i]][heatcool["j"][i]]*D.line_c[heatcool["i"][i]][heatcool["j"][i]]*nenh))
	if test<max_change*D.line_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=max_change*D.line_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	elif test>(1./max_change)*D.line_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=(1./max_change)*D.line_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	line_c_pre.append(test)	

	test=(heatcool["cool_ff"][i]/(D.brem_c_pre[heatcool["i"][i]][heatcool["j"][i]]*D.brem_c[heatcool["i"][i]][heatcool["j"][i]]*nenh))
	if test<max_change*D.brem_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=max_change*D.brem_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	elif test>(1./max_change)*D.brem_c_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=(1./max_change)*D.brem_c_pre[heatcool["i"][i]][heatcool["j"][i]]
	brem_c_pre.append(test)	

	test=(heatcool["heat_xray"][i]/(D.xray_h_pre[heatcool["i"][i]][heatcool["j"][i]]*D.xray_h[heatcool["i"][i]][heatcool["j"][i]]*nenh))
	if test<max_change*D.xray_h_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=max_change*D.xray_h_pre[heatcool["i"][i]][heatcool["j"][i]]
	elif test>(1./max_change)*D.xray_h_pre[heatcool["i"][i]][heatcool["j"][i]]:
		test=(1./max_change)*D.xray_h_pre[heatcool["i"][i]][heatcool["j"][i]]
	xray_h_pre.append(test)


fmt='%013.6e'

#This next line defines formats for the output variables. This is set in a dictionary
fmts={	'ir':'%03i',	
	'itheta':'%03i',	
	'rho':fmt,
	'comp_h_pre':fmt,
	'comp_c_pre':fmt,
	'xray_h_pre':fmt,
	'line_c_pre':fmt,
	'brem_c_pre':fmt,
	}	

titles=[]
titles=titles+["ir","itheta","rho"]
titles=titles+["comp_h_pre","comp_c_pre","xray_h_pre","brem_c_pre","line_c_pre"]	

col0=heatcool["i"]
col1=heatcool["j"]
col2=heatcool["rho"]
col3=comp_h_pre
col4=comp_c_pre
col5=xray_h_pre
col6=brem_c_pre
col7=line_c_pre

out=open("prefactors.dat",'w')

out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7],names=titles)
ascii.write(out_dat,out,formats=fmts)

