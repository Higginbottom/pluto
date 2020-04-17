#!/usr/bin/env python

import subprocess
import glob
from astropy.io import ascii
from astropy.table import Table
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import pyPLUTO as pp
import numpy as np



def pluto_input_file(tlim,data):
    output=open('pluto.ini','w')
    output.write("[Grid]\n")
    output.write("\n")
    if data["STRETCH_GRID"]:
        output.write("X1-grid 2 "+str(data["R_MIN"])+" "+str(data["N_R"])+" l+ "+str(data["R_MAX"])+" "+str(data["N_STRETCH"])+" s "+str(data["STRETCH_RMAX"])+"\n")
    else:
        output.write("X1-grid 2 "+str(data["R_MIN"])+" "+str(data["N_R"])+" l+ "+str(data["R_MAX"])+"\n")
    output.write("X2-grid 1    0.0    1    u    1.0\n")
    output.write("X3-grid 1    0.0    1    u    1.0\n")
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
    output.write("Solver         hll\n")
    output.write("\n")
    output.write("[Boundary]\n")
    output.write("\n")
    output.write("X1-beg        outflow\n")
    output.write("X1-end        outflow\n")
    output.write("X2-beg        outflow\n")
    output.write("X2-end        outflow\n")
    output.write("X3-beg        outflow\n")
    output.write("X3-end        outflow\n")
    output.write("\n")
    output.write("[Static Grid Output]\n")
    output.write("\n")
    if data["rad_force"]==1:
        output.write("uservar    4 dvdr M g_r T\n")
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
    output.write("R_0                         "+str(data["R_0"])+"\n")
    output.write("CENT_MASS                   "+str(data["CENT_MASS"])+"\n")
    output.write("T_ISO                       "+str(data["T_iso"])+"\n")
    output.write("V_0                         "+str(data["V_0"])+"\n")    
    output.write("M_rad                       "+str(data["M_rad"])+"\n")
    output.write("Lum_x                       "+str(data["Lum_x"])+"\n")
    output.write("Lum_UV                      "+str(data["Lum_UV"])+"\n")
    output.close()
    return

def pluto2py_1d(ifile):
    
    D=pp.pload(ifile)
    
    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')
    
    fname="%08d"%ifile+".pluto"
    out=open(fname,'w')
    
    
    ir=[]
    r=[]
    v_r=[]
    density=[]
    T=[]
    
    
    #The first element is a ghost cell
    ir.append(0)
    r.append(D.x1[0]*UNIT_LENGTH)
    v_r.append(D.vx1[0]*UNIT_VELOCITY)
    density.append(D.rho[0]*UNIT_DENSITY)
    T.append(D.T[0])
    
    for i in range(len(D.x1)):
        ir.append(i+1)
        r.append(D.x1[i]*UNIT_LENGTH)
        v_r.append(D.vx1[i]*UNIT_VELOCITY)
        density.append(D.rho[i]*UNIT_DENSITY)       
        T.append(D.T[i])
        
    #We need two extra ghost cells
    
    ir.append(i+1)
    r.append(D.x1[-1]*UNIT_LENGTH*1.1)
    v_r.append(D.vx1[-1]*UNIT_VELOCITY)
    density.append(D.rho[-1]*UNIT_DENSITY)
    T.append(D.T[-1])    
        
    ir.append(i+2)
    r.append(D.x1[-1]*UNIT_LENGTH*1.2)
    v_r.append(D.vx1[-1]*UNIT_VELOCITY)
    density.append(D.rho[-1]*UNIT_DENSITY)
    T.append(D.T[-1])
    
   
    fmt='%13.6e'

    #This next line defines formats for the output variables. This is set in a dictionary
    fmts={'ir':'%03i',    
        'r':fmt,
        'v_r':fmt,
        'density':fmt,
        'T':fmt}

    titles=["ir","r","v_r","density","T"]

    out_dat=Table([ir,r,v_r,density,T],names=titles)
    ascii.write(out_dat,out,formats=fmts)
    out.close()
     
    return(fname)
    

    
    
    
def python_input_file_stellar_wind(fname,data,cycles=2):
    output=open(fname+".pf",'w')
    output.write("System_type(star,binary,agn,previous)            star\n")
    output.write("\n")
    output.write("### Parameters for the Central Object\n")
    output.write("Central_object.mass(msol)                  "+str(data["CENT_MASS"])+"\n")
    output.write("Central_object.radius(cm)                  "+str(data["R_0"])+"\n")
    output.write("Central_object.radiation(yes,no)                  yes\n")
    output.write("Central_object.rad_type_to_make_wind(bb,models)                   bb\n")
    output.write("Central_object.temp                        "+str(data["T_star"])+"\n")    
    output.write("\n")
    output.write("### Parameters for the Disk (if there is one)\n")
    output.write("\n")
    output.write("Disk.type(none,flat,vertically.extended)       none\n")
    output.write("\n")
    output.write("### Parameters for BL or AGN\n")
    output.write("\n")
    output.write("Boundary_layer.radiation(yes,no)     no\n")
    output.write("\n")
    output.write("### Parameters descibing the various winds or coronae in the system\n")
    output.write("\n")
    output.write("Wind.radiation(yes,no) no\n")
    output.write("Wind.number_of_components  1\n")
    output.write("Wind.type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)  imported \n")
    output.write("Wind.coord_system(spherical,cylindrical,polar,cyl_var)  spherical\n")
    output.write("Wind.model2import              "+fname+"\n")
    output.write("\n")
    output.write("### Parameters associated with photon number, cycles,ionization and radiative transfer options\n")
    output.write("\n")
    output.write("Photons_per_cycle        "+str(data["NPHOT"])+"\n")
    output.write("Ionization_cycles        "+str(cycles)+"\n")
    output.write("Spectrum_cycles          0\n")
    output.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
    output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
    output.write("Atomic_data  data/standard80.dat\n")
    output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    reflect\n")
    output.write("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)   adiabatic\n")
    output.write("\n")
    output.write("### Parameters for Domain 0\n")
    output.write("\n")
    output.write("Wind.t.init                                40000\n")
    output.write("Wind.filling_factor(1=smooth,<1=clumped)   1\n")    
    output.write("Wind.radmax(cm)                            1e15\n")    
    output.write("\n")
    output.write("### Parameters for Reverberation Modeling (if needed)\n")
    output.write("\n")    
    output.write("Reverb.type(none,photon,wind,matom)   none\n")
    output.write("\n")    
    output.write("### Other parameters\n")
    output.write("\n")    
    output.write("Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)  logarithmic\n")
    output.write("Photon_sampling.nbands                     10\n")
    output.write("Photon_sampling.low_energy_limit(eV)       1.3333\n")
    output.write("Photon_sampling.high_energy_limit(eV)      200\n")
    
    output.close()
    return    


def driving_calc(ifile):


    max_accel_change=0.9

    D=pp.pload(ifile)

    flux=ascii.read("py_driving.dat")
    M=ascii.read("M_data.dat")
        

    print (flux.keys())

    # We need the definitions file - so we know the conversion factors.

    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')
    UNIT_ACCELERATION=UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH    
    
    gx_es=[]
    gy_es=[]
    gz_es=[]
    
    gx_bf=[]
    gy_bf=[]
    gz_bf=[]
    
    gx_line=[]
    gy_line=[]
    gz_line=[]
    
    gx=[]
    gy=[]
    gz=[]
    
    
    

    odd=0.0
    
    itest=0
    
    for i in range(len(flux["rho"])):
#        if (flux["rho"][i]/(D.rho[flux["i"][i]][flux["j"][i]]*UNIT_DENSITY))-1.>1e-6:
#            odd=odd+1
        

        #Electron scattering acceleration is taken directly from the python heatcool file
        gx_es.append(flux["rad_f_w"][i]/(flux["rho"][i]*flux["vol"][i]))
        gy_es.append(flux["rad_f_phi"][i]/(flux["rho"][i]*flux["vol"][i]))
        gz_es.append(flux["rad_f_z"][i]/(flux["rho"][i]*flux["vol"][i]))
        #BF scattering acceleration is taken directly from the python flux file
        gx_bf.append(flux["bf_f_w"][i]/(flux["rho"][i]*flux["vol"][i]))
        gy_bf.append(flux["bf_f_phi"][i]/(flux["rho"][i]*flux["vol"][i]))
        gz_bf.append(flux["bf_f_z"][i]/(flux["rho"][i]*flux["vol"][i]))
        #Line scattering acceleration is computed using the fluxes and the force multiplicayion factors
        gx_line.append((flux["F_vis_x"][i]*M["M_vis"][i]+flux["F_UV_x"][i]*M["M_uv"][i]+flux["F_Xray_x"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*flux["ne"][i]/flux["rho"][i]/c.c.cgs.value)
        gy_line.append((flux["F_vis_y"][i]*M["M_vis"][i]+flux["F_UV_y"][i]*M["M_uv"][i]+flux["F_Xray_y"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*flux["ne"][i]/flux["rho"][i]/c.c.cgs.value)
        gz_line.append((flux["F_vis_z"][i]*M["M_vis"][i]+flux["F_UV_z"][i]*M["M_uv"][i]+flux["F_Xray_z"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*flux["ne"][i]/flux["rho"][i]/c.c.cgs.value)
        #Add all the accelerations together
        gx.append(gx_es[-1]+gx_bf[-1]+gx_line[-1])
        gy.append(0.0)    
        gz.append(gz_es[-1]+gz_bf[-1]+gz_line[-1])


    
    fmt='%013.6e'

    #This next line defines formats for the output variables. This is set in a dictionary

    fmts2={'ir':'%03i',
        'rcent':fmt,
        'itheta':'%03i',
        'thetacent':fmt,    
        'gx':fmt,
        'gy':fmt,
        'gz':fmt,
        }

    titles=["ir","rcent","itheta","thetacent","rho"]
    titles=titles+["gx","gy","gz"]


    col0=flux["i"]
    col1=flux["rcen"]
    col2=flux["j"]
    col3=flux["thetacen"]
    col4=flux["rho"]

    out=open("py_accelerations.dat",'w')
    out_dat=Table([col0,col1,col2,col3,col4,gx,gy,gz],names=titles)
    ascii.write(out_dat,out,formats=fmts2)
    out.close()   
    

     
        
    return(odd)    



def get_units(fname='definitions.h'):
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
    inp.close()
    return(UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY)
    
def data_complete(data):    
    try:
        UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units()
    except:
        print("Unable to open definitions.h file - big problems")
        exit()
           
    
    #Rescale the grid
    
    data["R_MIN"]=data["R_MIN"]/UNIT_LENGTH
    data["R_MAX"]=data["R_MAX"]/UNIT_LENGTH
    data["STRETCH_RMAX"]=data["STRETCH_RMAX"]/UNIT_LENGTH
    
    return()
    