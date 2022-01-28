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
from scipy import interpolate


import pyPLUTO as pp
import numpy as np


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
    
    


def brem(freq,T_x=5.6e7,alpha=0.0):
    return freq**alpha*np.exp((-1.*c.h.cgs*freq/c.k_B.cgs/T_x).value)
    

    

def pluto_input_file(tlim,data):
    output=open('pluto.ini','w')
    output.write("[Grid]\n")
    output.write("\n")
    output.write("X1-grid 1 "+str(data["R_MIN"])+" "+str(data["N_R"])+" r "+str(data["R_MAX"])+" 1.07\n")
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
    output.write("Solver         hll\n")
    output.write("\n")
    output.write("[Boundary]\n")
    output.write("\n")
    output.write("X1-beg        reflective\n")
    output.write("X1-end        outflow\n")
    output.write("X2-beg        axisymmetric\n")
    output.write("X2-end        reflective\n")
    output.write("X3-beg        outflow\n")
    output.write("X3-end        outflow\n")
    output.write("\n")
    output.write("[Static Grid Output]\n")
    output.write("\n")
    if data["rad_force"]==1:
        output.write("uservar    20 XI T ch cc lc bc xh ch_pre cc_pre lc_pre bc_pre xh_pre ne nh gr gt gp dv_ds t M\n")
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
    output.write("MU                          %.2f"%(data["MU"])+"\n")
    output.write("RHO_0                       %4.2e"%(data["RHO_0"])+"\n")
    output.write("R_0                         %4.2e"%(data["R_0"])+"\n")
    output.write("RHO_ALPHA                   %.2f"%(data["RHO_ALPHA"])+"\n")
    output.write("CENT_MASS                   %4.2e"%(data["CENT_MASS"])+"\n")
    output.write("DISK_MDOT                   %4.2e"%(data["DISK_MDOT"])+"\n")
    output.write("T_ISO                       4e4  \n")
    output.write("L_star                      %4.2e"%(data["L_star"])+"\n")
    output.write("f_x                         %.2f"%(data["f_x"])+"\n") 
    output.write("f_uv                        %.2f"%(data["f_uv"])+"\n")
    output.write("T_x                         %4.2e"%(data["T_x"])+"\n")    
    output.write("KRAD                        %4.2e"%(data["k"])+"\n")
    output.write("ALPHARAD                    %4.2e"%(data["alpha"])+"\n")

#    output.write("L_x                         "+str(data["L_x"])+"\n")
#    output.write("T_x                         "+str(data["T_x"])+"\n")
#    output.write("DISK_TRUNC_RAD              "+str(data["DISK_TRUNC_RAD"])+"\n")
#    output.write("MU                          "+str(data["MU"])+"\n")
    output.close()
    return
    
    
    
    


def python_input_file(fname,data,cycles=2):
    output=open(fname+".pf",'w')
    output.write("System_type(star,binary,agn,previous)          "+data["system_type"]+"  \n")
    output.write("\n")
    output.write("### Parameters for the Central Object\n")
    output.write("Central_object.mass(msol)                  "+str(data["CENT_MASS"]/c.M_sun.cgs.value)+"\n")
    output.write("Central_object.radius(cm)                  "+str(data["CENT_RADIUS"])+"\n")
    output.write("\n")
    output.write("### Parameters for the Disk (if there is one)\n")
    output.write("\n")
    if data["disk_radiation"]=="yes":
        output.write("Disk.type(none,flat,vertically.extended)       flat\n")
        output.write("Disk.radiation(yes,no)      yes\n")
        output.write("Disk.temperature.profile(standard,readin,yso) standard\n") 
#        output.write("Disk.T_profile_file()  PSK_disk.dat\n")                
        output.write("Disk.rad_type_to_make_wind(bb,models) bb\n")        
        output.write("Disk.mdot(msol/yr) "+str(data["PY_DISK_MDOT"])+"\n")
        output.write("Disk.radmax(cm) "+str(data["DISK_TRUNC_RAD"])+"\n")
    else:
        output.write("Disk.radiation(yes,no)      no\n")        
    output.write("\n")
    output.write("### Parameters for BL or AGN\n")
    output.write("\n")
    if data["boundary_layer"]=="no":
        output.write("Boundary_layer.radiation(yes,no)                   no\n")
    else:
        output.write("Boundary_layer.radiation(yes,no)                   yes\n")        
        output.write("Boundary_layer.rad_type_to_make_wind(bb,models,power) bb\n")
        output.write("Boundary_layer.temp(K) "+str(data["T_BL"])+"\n")
        output.write("Boundary_layer.luminosity(ergs/s)  "+str(data["L_BL"])+"\n")        
    if data["cent_spectype"]=="brem":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   brems\n")
        output.write("AGN.bremsstrahlung_temp(K) "+str(data["T_x"])+"\n")
        output.write("Central_object.luminosity(ergs/s) "+str(data["L_2_10"])+"\n")
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")
    elif data["cent_spectype"]=="bb":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   bb\n")
        output.write("Central_object.blackbody_temp(K)                        "+str(data["T_star"])+"\n") 
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")   
    elif data["cent_spectype"]=="models":
        output.write("Central_object.radiation(yes,no)     yes\n")
        output.write("Central_object.rad_type_to_make_wind(bb,models,power,cloudy,brems)   models\n")
        output.write("Input_spectra.model_file          model.ls\n")
        output.write("Central_object.luminosity(ergs/s) "+str(data["L_2_10"])+"\n")
        output.write("Central_object.geometry_for_source(sphere,lamp_post) sphere\n")
    elif data["cent_spectype"]=="none":
        output.write("Central_object.radiation(yes,no)     no\n")
    output.write("\n")
    output.write("### Parameters descibing the various winds or coronae in the system\n")
    output.write("\n")
    if data["wind_radiation"]=="yes":
        output.write("Wind.radiation(yes,no) yes\n")
    else:
        output.write("Wind.radiation(yes,no) no\n")        
    output.write("Wind.number_of_components  1\n")
    output.write("Wind.type(SV,star,hydro,corona,kwd,homologous,yso,shell,imported)  imported \n")
    output.write("Wind.coord_system(spherical,cylindrical,polar,cyl_var)  polar\n")
    output.write("Wind.dim.in.x_or_r.direction               30\n")
    output.write("Wind.dim.in.z_or_theta.direction           30\n")
    output.write("\n")
    output.write("### Parameters associated with photon number, cycles,ionization and radiative transfer options\n")
    output.write("\n")
    output.write("Photons_per_cycle        "+str(data["NPHOT"])+"\n")
    output.write("Ionization_cycles        "+str(cycles)+"\n")
    output.write("Spectrum_cycles          0\n")
    output.write("Wind.ionization(on.the.spot,ML93,LTE_tr,LTE_te,fixed,matrix_bb,matrix_pow)  matrix_pow\n")
    if data["line_trans"]=="macro":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   macro_atoms_escape_prob\n")
        output.write("Atomic_data  data/h10_hetop_standard80.dat\n")
        output.write("Matom_transition_mode(mc_jumps,matrix) matrix\n")
    elif data["line_trans"]=="simple":
        output.write("Line_transfer(pure_abs,pure_scat,sing_scat,escape_prob,thermal_trapping,macro_atoms,macro_atoms_thermal_trapping)   escape_prob\n")
        output.write("Atomic_data  data/standard80.dat\n")
    output.write("Surface.reflection.or.absorption(reflect,absorb,thermalized.rerad)    thermalized.rerad\n")
    output.write("Wind_heating.extra_processes(none,adiabatic,nonthermal,both)   none\n")
    output.write("\n")
    output.write("### Parameters for Domain 0\n")
    output.write("\n")
    output.write("Wind.model2import "+fname+".pluto\n")
#    output.write("Hydro.thetamax(degrees:negative_means_no_maximum)  -1\n")
    output.write("Wind.t.init                                40000\n")
    output.write("Wind.filling_factor(1=smooth,<1=clumped)   1\n")
    output.write("\n")
    output.write("### Parameters for Reverberation Modeling (if needed)\n")
    output.write("\n")    
    output.write("Reverb.type(none,photon,wind,matom)   none\n")
    output.write("\n")    
    output.write("### Other parameters\n")
    output.write("\n")    
    output.write("Photon_sampling.approach(T_star,cv,yso,AGN,min_max_freq,user_bands,cloudy_test,wide,logarithmic)  logarithmic\n")
    output.write("Photon_sampling.nbands                     10\n")
    output.write("Photon_sampling.low_energy_limit(eV)       0.13333\n")
    output.write("Photon_sampling.high_energy_limit(eV)      100000\n")
    
    output.close()
    return
    
    

def pluto2py_rtheta(ifile):



    D=pp.pload(ifile)
    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')



    pluto_r_inner=D.x1r*UNIT_LENGTH #Python expects a model file to have the inner radius of each shell these are the pluto coordinated in the inner edge
    pluto_theta_inner=D.x2r #Theta coordinates of inner edges of cell


    pluto_r_c=D.x1*UNIT_LENGTH #Central point of each pluto cell - where the velocity is defined
    pluto_theta_c=D.x2  #Central theta of each pluto cell - where the velocity is defined

    pluto_vr_c=D.vx1*UNIT_VELOCITY #The velocity in each pluto cell - defined at cell centre
    pluto_vt_c=D.vx2*UNIT_VELOCITY #The velocity in each pluto cell - defined at cell centre
    pluto_vy_c=D.vx3*UNIT_VELOCITY #The y component is just equal to the phi component....

    pluto_density=D.rho*UNIT_DENSITY
    pluto_temperature=D.T


    pluto_vx_c=np.zeros(np.shape(pluto_vr_c))
    pluto_vz_c=np.zeros(np.shape(pluto_vr_c))


    python_r=[]
    python_theta=[]

    for i in range(len(pluto_r_inner)):
        python_r.append(pluto_r_inner[i])
    python_r.append(pluto_r_inner[-1]+(pluto_r_inner[-1]-pluto_r_inner[-2]))
    
    for i in range(len(pluto_theta_inner)):
        python_theta.append(pluto_theta_inner[i])
    python_theta.append(pluto_theta_inner[-1]+(pluto_theta_inner[-1]-pluto_theta_inner[-2]))



    for i in range(len(pluto_r_c)):
        for j in range(len(pluto_theta_c)):
                pluto_vx_c[i][j]=pluto_vr_c[i][j]*np.sin(pluto_theta_c[j])+pluto_vt_c[i][j]*np.cos(pluto_theta_c[j])
                pluto_vz_c[i][j]=pluto_vr_c[i][j]*np.cos(pluto_theta_c[j])-pluto_vt_c[i][j]*np.sin(pluto_theta_c[j])



    #First, we interpolate in the r direction

    vx_temp=np.zeros([len(python_r),len(pluto_theta_c)])
    vy_temp=np.zeros([len(python_r),len(pluto_theta_c)])
    vz_temp=np.zeros([len(python_r),len(pluto_theta_c)])

    for i in range(len(pluto_theta_c)):
        vx=interpolate.interp1d(pluto_r_c,pluto_vx_c[:,i],fill_value='extrapolate')
        vy=interpolate.interp1d(pluto_r_c,pluto_vy_c[:,i],fill_value='extrapolate')
        vz=interpolate.interp1d(pluto_r_c,pluto_vz_c[:,i],fill_value='extrapolate')
        vx_temp[:,i]=vx(python_r)
        vy_temp[:,i]=vy(python_r)
        vz_temp[:,i]=vz(python_r)
    

    #And now we interpolate in the theta direction

    python_vx=np.zeros([len(python_r),len(python_theta)])
    python_vy=np.zeros([len(python_r),len(python_theta)])
    python_vz=np.zeros([len(python_r),len(python_theta)])
 
    for i in range(len(python_r)):
        vx=interpolate.interp1d(pluto_theta_c,vx_temp[i],fill_value='extrapolate')
        vy=interpolate.interp1d(pluto_theta_c,vy_temp[i],fill_value='extrapolate')
        vz=interpolate.interp1d(pluto_theta_c,vz_temp[i],fill_value='extrapolate')
        python_vx[i]=vx(python_theta)
        python_vy[i]=vy(python_theta)
        python_vz[i]=vz(python_theta) 
    
    #Finally deal with cell centred values, rho and T   

    python_density=np.zeros([len(python_r),len(python_theta)])
    python_T_e=np.zeros([len(python_r),len(python_theta)])

    python_density[0:-2,0:-2]=pluto_density
    python_T_e[0:-2,0:-2]=pluto_temperature



    #Now we need to turn these arrays into linear vectors for output

    ir=[]
    r=[]
    itheta=[]
    theta=[]
    inwind=[]
    v_x=[]
    v_y=[]
    v_z=[]
    density=[]
    T=[]



    for i in range(len(python_r)):
        for j in range(len(python_theta)):
            ir.append(i)
            r.append(python_r[i])
            itheta.append(j)
            theta.append(np.degrees(python_theta[j]))
            if python_density[i][j]==0.0:
                inwind.append(-1)   #outside wind
            else:
                inwind.append(0)    #Fully in wind
            v_x.append(python_vx[i][j])
            v_y.append(python_vy[i][j])
            v_z.append(python_vz[i][j])
            density.append(python_density[i][j])
            T.append(python_T_e[i][j])
        



    titles=["ir","itheta","inwind","r","theta","v_x","v_y","v_z","density","T"]

    #This next line defines formats for the output variables. This is set in a dictionary
    fmt='%13.6e'
    fmts={'ir':'%03i',
        'itheta':'%03i',   
        'inwind':'%01i',         
        'r':fmt,
        'theta':fmt,
        'v_x':fmt,
        'v_y':fmt,
        'v_z':fmt,
        'density':fmt,
        'T':fmt}



    fname="%08d"%ifile+".pluto"
    out=open(fname,'w')

    out_dat=Table([ir,itheta,inwind,r,theta,v_x,v_y,v_z,density,T],names=titles)
    ascii.write(out_dat,out,formats=fmts)
    out.close()
    return(fname)


        
    
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
    out.write("#                r, theta phi for spherical polars\n")
    out.write("#                         x,y,z        for carteisan\n")
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
    fmts={    'ir':'%03i',    
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
    max_change=0.1
    max_accel_change=0.9

    heatcool=ascii.read("%08d"%(ifile)+"_py_heatcool.dat")
    D=pp.pload(ifile)

    if radforce:
        py_driving=ascii.read("%08d"%(ifile)+"_py_driving.dat")
        M=ascii.read("%08d"%(ifile)+"_M_data.dat")
        

    # We need the definitions file - so we know the conversion factors.

    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')
    UNIT_ACCELERATION=UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH    

    comp_h_pre=[]
    comp_c_pre=[]
    xray_h_pre=[]
    brem_c_pre=[]
    line_c_pre=[]
    
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
    
        ideal_prefactor=(heatcool["heat_xray"][i]/(D.xh[heatcool["i"][i]][heatcool["j"][i]]*nhnh))
        change=ideal_prefactor/D.xh_pre[heatcool["i"][i]][heatcool["j"][i]]
        if change<max_change:
            change=max_change
        elif change>(1./max_change):
            change=(1./max_change)
        xray_h_pre.append(change*D.xh_pre[heatcool["i"][i]][heatcool["j"][i]])
        
        if radforce:
            #Electron scattering acceleration is taken directly from the python heatcool file
            gx_es.append(py_driving["es_f_x"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gy_es.append(py_driving["es_f_y"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gz_es.append(py_driving["es_f_z"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            #BF scattering acceleration is taken directly from the python py_driving file
            gx_bf.append(py_driving["bf_f_x"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gy_bf.append(py_driving["bf_f_y"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gz_bf.append(py_driving["bf_f_z"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            #line scattering acceleration is taken directly from the python py_driving file
            gx_line.append(py_driving["ka_f_x"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gy_line.append(py_driving["ka_f_y"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gz_line.append(py_driving["ka_f_z"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            #Line scattering acceleration is computed using the fluxes and the force multiplicayion factors
#  gx_line.append((py_driving["F_vis_x"][i]*M["M_vis"][i]+py_driving["F_UV_x"][i]*M["M_uv"][i]+py_driving["F_Xray_x"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)
#            gy_line.append((py_driving["F_vis_y"][i]*M["M_vis"][i]+py_driving["F_UV_y"][i]*M["M_uv"][i]+py_driving["F_Xray_y"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)
#            gz_line.append((py_driving["F_vis_z"][i]*M["M_vis"][i]+py_driving["F_UV_z"][i]*M["M_uv"][i]+py_driving["F_Xray_z"][i]*M["M_xray"][i])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)

            #Add all the accelerations together
            gx.append(gx_es[-1]+gx_bf[-1]+gx_line[-1])
            gy.append(0.0)    
            gz.append(gz_es[-1]+gz_bf[-1]+gz_line[-1])
#            gx.append(heatcool["rad_f_w"][i]/heatcool["rho"][i]/heatcool["vol"][i])    
#            gy.append(0.0)    
#            gz.append(heatcool["rad_f_z"][i]/heatcool["rho"][i]/heatcool["vol"][i])    
        else:        
            gx.append(0.0)
            gy.append(0.0)
            gz.append(0.0)
            
    
    fmt='%013.6e'

    #This next line defines formats for the output variables. This is set in a dictionary
    fmts={    'ir':'%03i',
        'rcent':fmt,
        'itheta':'%03i',
        'thetacent':fmt,    
        'rho':fmt,
        'comp_h_pre':fmt,
        'comp_c_pre':fmt,
        'xray_h_pre':fmt,
        'line_c_pre':fmt,
        'brem_c_pre':fmt,
        'gx':fmt,
        'gy':fmt,
        'gz':fmt,
        }  
          
    titles=[]
    titles=titles+["ir","rcent","itheta","thetacent","rho"]
    titles=titles+["comp_h_pre","comp_c_pre","xray_h_pre","brem_c_pre","line_c_pre"]
    titles=titles+["gx","gy","gz"]    
    
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
    col10=gx
    col11=gy
    col12=gz

    out=open("prefactors.dat",'w')

    out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12],names=titles)
    ascii.write(out_dat,out,formats=fmts)
    out.close()
    
    if radforce:
        fmts2={'ir':'%03i',
            'rcent':fmt,
            'itheta':'%03i',
            'thetacent':fmt,    
            'gx_es':fmt,
            'gy_es':fmt,
            'gz_es':fmt,
            'gx_bf':fmt,
            'gx_bf':fmt,
            'gx_bf':fmt,
            'gx_line':fmt,
            'gx_line':fmt,
            'gx_line':fmt,
            }
    
        titles=[]
        titles=titles+["ir","rcent","itheta","thetacent"]
        titles=titles+["gx_es","gy_es","gz_es"]
        titles=titles+["gx_bf","gy_bf","gz_bf"]
        titles=titles+["gx_line","gy_line","gz_line"]
    
        out=open("accelerations.dat",'w')
        out_dat=Table([col0,col1,col2,col3,gx_es,gy_es,gz_es,gx_bf,gy_bf,gz_bf,gx_line,gy_line,gz_line],names=titles)
        ascii.write(out_dat,out,formats=fmts2)
        out.close()    
    
    
    return(odd)
        
def pre_calc_k_alpha(ifile,radforce=0):
    max_change=0.9
    max_accel_change=0.9
    
    k=0.2
    alpha=0.6
    
    heatcool=ascii.read("%08d"%(ifile)+"_py_heatcool.dat")
    D=pp.pload(ifile)

    if radforce:
        py_driving=ascii.read("%08d"%(ifile)+"_py_driving.dat")
        py_pcon=ascii.read("%08d"%(ifile)+"_py_pcon_data.dat",data_start=1)
        t_vis=py_pcon['col7']
        t_UV=py_pcon['col8']
        t_Xray=py_pcon['col9']
        

    # We need the definitions file - so we know the conversion factors.

    UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units('definitions.h')
    UNIT_ACCELERATION=UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH    

    comp_h_pre=[]
    comp_c_pre=[]
    xray_h_pre=[]
    brem_c_pre=[]
    line_c_pre=[]
    
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
    
    M_vis=[]
    M_uv=[]
    M_xray=[]
    M_max=4400.
    
    eta_max=(M_max/(k*(1.-alpha)))**(1./alpha)
    for i in range(len(heatcool["rho"])):
        
        comp_h_pre.append(0.0)
        comp_c_pre.append(0.0)
        line_c_pre.append(0.0)
        brem_c_pre.append(0.0)
        xray_h_pre.append(0.0)
        

        
        
        if radforce:

#            M_v=(np.min([k*t_vis[i]**(-1.0*alpha),M_max]))
#            M_u=(np.min([k*t_UV[i]**(-1.0*alpha),M_max]))
#            M_x=(np.min([k*t_Xray[i]**(-1.0*alpha),M_max]))
            
            
            tau_max=t_vis[i]*eta_max 
            
            if (t_vis[i]>0.0):
                M_v=k*t_vis[i]**(-1.0*alpha)*((1.+tau_max)**(1.-alpha)-1.)/(tau_max**(1.-alpha))
            else:
                M_v=0.0
            tau_max=t_UV[i]*eta_max
            if (t_UV[i]>0.0):
                M_u=k*t_UV[i]**(-1.0*alpha)*((1.+tau_max)**(1.-alpha)-1.)/(tau_max**(1.-alpha))
            else:
                M_u=0.0
            tau_max=t_Xray[i]*eta_max
            if (t_Xray[i]>0.0):
                M_x=k*t_Xray[i]**(-1.0*alpha)*((1.+tau_max)**(1.-alpha)-1.)/(tau_max**(1.-alpha))
            else:
                M_x=0.0

            M_vis.append(M_v+1.)
            M_uv.append(M_u+1.)
            M_xray.append(M_x+1.)
            
            
            
            #Electron scattering acceleration is taken directly from the python heatcool file
            gx_es.append(py_driving["es_f_x"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gy_es.append(py_driving["es_f_y"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gz_es.append(py_driving["es_f_z"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            #BF scattering acceleration is taken directly from the python py_driving file
            gx_bf.append(py_driving["bf_f_x"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gy_bf.append(py_driving["bf_f_y"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            gz_bf.append(py_driving["bf_f_z"][i]/(py_driving["rho"][i]*py_driving["vol"][i]))
            #Line scattering acceleration is computed using the fluxes and the force multiplicayion factors
            gx_line.append((py_driving["F_vis_x2"][i]*M_vis[-1]+py_driving["F_UV_x2"][i]*M_uv[-1]+py_driving["F_Xray_x2"][i]*M_xray[-1])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)
            gy_line.append((py_driving["F_vis_y2"][i]*M_vis[-1]+py_driving["F_UV_y2"][i]*M_uv[-1]+py_driving["F_Xray_y2"][i]*M_xray[-1])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)
            gz_line.append((py_driving["F_vis_z2"][i]*M_vis[-1]+py_driving["F_UV_z2"][i]*M_uv[-1]+py_driving["F_Xray_z2"][i]*M_xray[-1])*c.sigma_T.cgs.value*py_driving["ne"][i]/py_driving["rho"][i]/c.c.cgs.value)

            #Add all the accelerations together
            gx.append(gx_bf[-1]+gx_line[-1])
            gy.append(0.0)    
            gz.append(gz_bf[-1]+gz_line[-1])
#            gx.append(heatcool["rad_f_w"][i]/heatcool["rho"][i]/heatcool["vol"][i])    
#            gy.append(0.0)    
#            gz.append(heatcool["rad_f_z"][i]/heatcool["rho"][i]/heatcool["vol"][i])    
        else:        
            gx.append(0.0)
            gy.append(0.0)
            gz.append(0.0)
            
    
    fmt='%013.6e'

    #This next line defines formats for the output variables. This is set in a dictionary
    fmts={    'ir':'%03i',
        'rcent':fmt,
        'itheta':'%03i',
        'thetacent':fmt,    
        'rho':fmt,
        'comp_h_pre':fmt,
        'comp_c_pre':fmt,
        'xray_h_pre':fmt,
        'line_c_pre':fmt,
        'brem_c_pre':fmt,
        'gx':fmt,
        'gy':fmt,
        'gz':fmt,
        }  
          
    titles=[]
    titles=titles+["ir","rcent","itheta","thetacent","rho"]
    titles=titles+["comp_h_pre","comp_c_pre","xray_h_pre","brem_c_pre","line_c_pre"]
    titles=titles+["gx","gy","gz"]    
    
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
    col10=gx
    col11=gy
    col12=gz

    out=open("prefactors.dat",'w')

    out_dat=Table([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12],names=titles)
    ascii.write(out_dat,out,formats=fmts)
    out.close()
    
    
    
    
    if radforce:
        fmts2={'ir':'%03i',
            'rcent':fmt,
            'itheta':'%03i',
            'thetacent':fmt,  
            'rho':fmt,     
            'gx':fmt,
            'gy':fmt,
            'gz':fmt,
            }
    
        titles=[]
        titles=titles+["ir","rcent","itheta","thetacent","rho"]

        titles=titles+["gx","gy","gz"]
    
        out=open("py_accelerations.dat",'w')
        out_dat=Table([col0,col1,col2,col3,col4,gx,gy,gz],names=titles)
        ascii.write(out_dat,out,formats=fmts2)
        out.close() 
   
   
   
        
        fmts3={'ir':'%03i',
            'rcent':fmt,
            'itheta':'%03i',
            'thetacent':fmt, 
            'M_vis':fmt,
            'M_uv':fmt,
            'M_xray':fmt
            }
        
        titles=[]
        titles=titles+["ir","rcent","itheta","thetacent"]
        titles=titles+["M_vis","M_uv","M_xray"]
        out=open("M_k_alpha_data.dat",'w')
        out_dat=Table([col0,col1,col2,col3,M_vis,M_uv,M_xray],names=titles)
        ascii.write(out_dat,out,formats=fmts3)
        out.close()                    
    
    
    return(odd)    
    
def data_complete(data):    
    try:
        UNIT_DENSITY,UNIT_LENGTH,UNIT_VELOCITY=get_units()
    except:
        print("Unable to open definitions.h file - big problems")
        exit()
    data["efficiency"]=0.083        #The efficiency of conversion of mass to lumonisity at the central source - used to set the disk mdot
    data["DISK_MDOT"]=(data["L_x"]/c.c.cgs/c.c.cgs/data["efficiency"]).value   #THe disk massloss rate is only used to set the initial temperature
    data["PY_DISK_MDOT"]=data["DISK_MDOT"]*365.25*60*60*24/c.M_sun.cgs.value #The disk massloss rate for python
    
    #Now we work out the matching luminosity for the python simulation.

    nu1=((13.6*u.eV).to(u.erg)/c.h.cgs).value
    nu2=((2000*u.eV).to(u.erg)/c.h.cgs).value
    nu3=((10000*u.eV).to(u.erg)/c.h.cgs).value
    
    if data["spectype"]=='brem':    
        numax=(data["T_x"]*u.K*c.k_B.cgs/c.h.cgs).value*100.
        L_x_test=quad(brem,nu1,np.max([nu3,numax]),args=(data["T_x"],data["BREM_ALPHA"]))
        const=data["L_x"]/L_x_test[0]
        data["L_2_10"]=const*quad(brem,nu2,nu3,args=(data["T_x"],data["BREM_ALPHA"]))[0]    
    elif data["spectype"]=="models":
        try:
            inp=open(data["model"])
            mod_file=inp.readline().split()[0]
            inp.close()
        except:
            print("Cannot open model file")
            exit()
        mod_lambda=[]
        mod_llambda=[]
        print("opening model file",mod_file)
        inp=open(mod_file)
        for line in inp.readlines():
            data_split=line.split()
            mod_lambda.append(float(data_split[0]))
            mod_llambda.append(float(data_split[1]))
        inp.close()
        mod_lambda=np.array(mod_lambda)
        mod_llambda=np.array(mod_llambda)
        mod_freq=c.c.cgs.value/(mod_lambda[::-1]/1e8)
        mod_lnu=(mod_llambda/(c.c.cgs.value*1e8/mod_lambda/mod_lambda))[::-1]
        
        xvals=np.linspace(nu1,mod_freq[-1],1e6)
        yvals=np.interp(xvals,mod_freq,mod_lnu)
        L_x=np.trapz(yvals,xvals)
        
        xvals=np.linspace(nu2,nu3,1e6)
        yvals=data["L_x"]/L_x*np.interp(xvals,mod_freq,mod_lnu)
        data["L_2_10"]=np.trapz(yvals,xvals)
            
    
    #Rescale the grid
    
    data["R_MIN"]=data["R_MIN"]/UNIT_LENGTH
    data["R_MAX"]=data["R_MAX"]/UNIT_LENGTH
    
    return()

def get_status():
    cmdline="tail -1 dbl.out"   
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) 
    dbl_file_1=int(proc.stdout.read().split()[0])
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE)
    last_dbl_time=float(proc.stdout.read().split()[1])
    print("Last file listed in dbl file:",dbl_file_1)
    cmdline="ls  *.dbl | tail -1"
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) 
    dbl_file_2=int(proc.stdout.read().split(b'.')[1])
    print("last dbl file in directory:  ",dbl_file_2)    
    cmdline="grep tstop pluto.ini"
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) 
    dbl_time_requested=float(proc.stdout.read().split()[1])
    print ("Last dbl file was for time:  ",last_dbl_time)
    print ("Last dbl file requested time:",dbl_time_requested)
    cmdline="tail -50 input.sig | grep Finished | grep cycle | tail -1"
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) 
    test=proc.stdout.read().split()
    py_last_completed=int(test[8])
    cmdline="tail -50 input.sig | grep Starting | grep cycle | tail -1"
    proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) 
    test=proc.stdout.read().split()
    py_last_requested=int(test[10])
    return dbl_file_1,dbl_file_2,last_dbl_time,dbl_time_requested,py_last_completed,py_last_requested 
    
    
    
       
    
def python_input_gen(ifile,py_cycles,data):
    pluto2py(ifile)   #We now make a python input file
    root="%08d"%(ifile)
    python_input_file(root+".pluto",data,py_cycles)  #This generate a python parameter file
    cmdline="cp "+root+".pluto"+".pf input.pf"   #Copy the python file to a generaic name so windsave files persist
    print((cmdline+"\n"))
    subprocess.check_call(cmdline,shell=True)
    return
    
def pluto_input_gen(ifile,data):
    root="%08d"%(ifile)
    cmdline="cp py_heatcool.dat "+root+"_py_heatcool.dat"  
    subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
    if data["rad_force"]:
        print ("Running CAK\n")
        cmdline=" ./cak > cak_output" 
        print (cmdline)         
        subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
#    now make a prefactors file
    print ("Making a prefactor file using "+str(ifile)+" dbl file")
    pre_calc(ifile,data["rad_force"])    
    cmdline="cp prefactors.dat "+root+"_prefactors.dat"  
    subprocess.check_call(cmdline,shell=True)
    cmdline="cp input.wind_save "+root+"_input.wind_save"  
    subprocess.check_call(cmdline,shell=True)
    if ifile>2:
        cmdline="rm "+"%08d"%(ifile-2)+"_input.wind_save"
        print (cmdline+"\n")
        try:
            subprocess.check_call(cmdline,shell=True)
        except:
            print("Could not delete\n")
    
    
    
def loop(t0,dt,istart,py_cycles,data,flag):
    
    ifile=istart #If this is a restart - we need to populate ifile
    
    root="%08d"%(istart)
    out=open("pluto_py_logfile",'w')
    out.write("Starting run"+"\n")
    #out.write("zeus_ver="+zeus_ver+"\n")
    for i in range(istart,10000):  
        out.write("STARTING CYCLE "+str(i)+"\n")
        print(("STARTING CYCLE "+str(i)+"\n"))
        if flag==0: #If flag is not set - then we run pluto first 
            if i==0: #We are on our first cycle
                time=t0
            else: #We have already done at least one pluto run - we need to increase the time
                 cmdline="tail -1 dbl.out"   
                 proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE)
                 last_dbl_time=float(proc.stdout.read().split()[1])
                 time=last_dbl_time+dt
            pluto_input_file(time,data)
            out.write("Running for time="+str(time)+"\n")
            print("Running for time="+str(time)+"\n")
            
            if i==0:   #This is the first step - 
                out.write("Creating first zeus_file"+"\n")
                cmdline="mpirun -n "+str(data["nproc_pl"])+" ./pluto >"+"%08d"%i+"_pluto_log"
            else:
                out.write("generating restart zeus run \n")      #This should be the name of the restart file
                cmdline="mpirun -n "+str(data["nproc_pl"])+" ./pluto -restart "+str(ifile)+" > "+"%08d"%i+"_pluto_log"
            print("Executing pluto with command line "+cmdline)
            out.write("Executing pluto with command line "+cmdline+"\n")
            return_number=subprocess.call(cmdline,shell=True)    #Call pluto
            if return_number==0:
                out.write("Finished pluto run successfully"+"\n")
            else:    
                out.write("Pluto run crashed"+"\n")
                exit(0)
            cmdline="tail -1 dbl.out"   
            out.write(cmdline+"\n")
            proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) #This mess gets the last dblfile    
            ifile=int(proc.stdout.read().split()[0])
            pluto2py(ifile)   #We now make a python input file
            root="%08d"%(ifile)
            py_cycles=py_cycles+2
            python_input_file(root+".pluto",data,py_cycles)  #This generate a python parameter file
            cmdline="cp "+root+".pluto"+".pf input.pf"   #Copy the python file to a generaic name so windsave files persist
            out.write(cmdline+"\n")
            print((cmdline+"\n"))
            subprocess.check_call(cmdline,shell=True)
            
        
        if py_cycles==3: #This is the first time through - so no restart""
            cmdline="mpirun -n "+str(data["nproc_py"])+" "+data["python_ver"]+" -z  input.pf > "+root+".python_log"  #We now run python
        else:
            cmdline="mpirun -n "+str(data["nproc_py"])+" "+data["python_ver"]+" -z -r  input.pf > "+root+".python_log"  #We now run python
        out.write("Running python"+"\n") 
        print(("Running python"+"\n"))     
        out.write(cmdline+"\n")
        print((cmdline+"\n"))
        subprocess.check_call(cmdline,shell=True)   #Well, here is the actual call
        cmdline="cp py_heatcool.dat "+root+"_py_heatcool.dat"  
        out.write(cmdline+"\n")
        subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
        flag=0 #reset the flag that may have been set on entry to run python first
        if data["rad_force"]:
            out.write ("Running CAK\n")
            cmdline="mpirun -n "+str(data["nproc_cak"])+" ./cak > cak_output" 
            out.write (cmdline)
             
            subprocess.check_call(cmdline,shell=True)   #And finally we take a copy of the python heatcool file for later investigation.
    #    now make a prefactors file
        out.write ("Making a prefactor file using "+str(ifile)+" dbl file")
        pre_calc(ifile,data["rad_force"])    
        cmdline="cp prefactors.dat "+root+"_prefactors.dat"  
        out.write(cmdline+"\n")
        subprocess.check_call(cmdline,shell=True)
        if data["rad_force"]:
            cmdline="cp accelerations.dat "+root+"_accelerations.dat"  
            out.write(cmdline+"\n")
            subprocess.check_call(cmdline,shell=True)
        cmdline="cp input.wind_save "+root+"_input.wind_save"  
        out.write(cmdline+"\n")
        subprocess.check_call(cmdline,shell=True)
        if ifile>2:
            cmdline="rm "+"%08d"%(ifile-2)+"_input.wind_save"
            out.write(cmdline+"\n")
            try:
                subprocess.check_call(cmdline,shell=True)
            except:
                out.write("Could not delete\n")
                
                
            
        out.write("FINISHED CYCLE"+"\n")
    out.close()
    
