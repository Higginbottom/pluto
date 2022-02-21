#include "pluto.h"

double heatcool2();
double sqsqxi;

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv,ii,iangle;  
#if COOLING != NO 
double ***comp_h_pre, ***comp_c_pre, ***line_c_pre, ***brem_c_pre, ***xray_h_pre;
#endif

double ***T_out, ***xi_out, ***ne_out, ***nH_out,krad,alpharad,krad_UV,alpha_UV,t_UV,M_UV;
double M_UV_array[MPOINTS];


#if (BODY_FORCE & VECTOR)
double ***gr_out, ***gt_out, ***gp_out, g[3], Vc[NVAR];
double ***gx_pre, ***gy_pre, ***gz_pre;
#endif

double ***dv_ds_out;
double ***t_out,***M_out1,***M_out2,M_max,eta_max,tau_max,M_temp1,M_temp2,fmax;

double rho,p,T,xi,ne,nH,n,mu,lx,tx,r,sigma_e,v_th;
double *x1   = grid->x[IDIR],  *x2  = grid->x[JDIR],  *x3  = grid->x[KDIR];


#if COOLING != NO 
comp_h_pre  = GetUserVar("ch_pre");
comp_c_pre  = GetUserVar("cc_pre");
xray_h_pre  = GetUserVar("xh_pre");
line_c_pre  = GetUserVar("lc_pre");
brem_c_pre  = GetUserVar("bc_pre");
#endif

ne_out	= GetUserVar("ne");
nH_out	= GetUserVar("nh");
T_out     = GetUserVar("T");
xi_out     = GetUserVar("XI");

#if (BODY_FORCE & VECTOR)
gr_out     = GetUserVar("gr");
gt_out     = GetUserVar("gt");
gp_out     = GetUserVar("gp");
#endif


dv_ds_out=GetUserVar("dv_ds");
t_out=GetUserVar("t");
M_out1=GetUserVar("M_max1");
M_out2=GetUserVar("M_max2");





mu=g_inputParam[MU];
lx=g_inputParam[L_star]*g_inputParam[f_x];  //Xray luminosiy
tx=g_inputParam[T_x];  //Xray tenperature 

  
sigma_e=CONST_sigmaT/CONST_amu/1.18;
  
  
  DOM_LOOP(k,j,i){
  	  r=grid->x[IDIR][i]*UNIT_LENGTH;  //The radius - in real units
      rho = d->Vc[RHO][k][j][i]*UNIT_DENSITY;  //density of the current cell in physical units	 
      
      #if EOS==ISOTHERMAL
      T = T_out[k][j][i]=g_inputParam[T_ISO];
      #else
      T   = T_out[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*mu;    //Compute initial temperature in Kelvin
      #endif
	  nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming stellar abundances	  
  	  xi_out[k][j][i]=xi=lx/nH/r/r;     //ionization parameter	  
	  
	  n=rho/(mu*CONST_mp);    //particle density	  
      nH_out[k][j][i]=nH; 
	  ne=1.21*nH;             //electron number density assuming full ionization
	  
	  
	  
      v_th=  pow ((2. * CONST_kB * T / CONST_mp), 0.5);     //We need the thermal velocity for hydrogen
	  
#if COOLING==BLONDIN
	 
	 heatcool2(xi,T,i,j,k,ne,nH);
     ne_out[k][j][i]=ne;	  
	comp_h_pre[k][j][i]=d->comp_h_pre[k][j][i];
	comp_c_pre[k][j][i]=d->comp_c_pre[k][j][i];
	xray_h_pre[k][j][i]=d->xray_h_pre[k][j][i];
	line_c_pre[k][j][i]=d->line_c_pre[k][j][i];
	brem_c_pre[k][j][i]=d->brem_c_pre[k][j][i];
#endif
/* If we have a vector pody force we need to popkluate and output the relevent arrays.
	Note that there is a complecity here - the forces are in sphperical polar, but the prefactors are in cartesian
	This is to allow us to deal with non radial forces. */	
	
	
	#if (BODY_FORCE & VECTOR)
	NVAR_LOOP(nv)  Vc[nv]=d->Vc[nv][k][j][i];	
	
	BodyForceVector(Vc, g, x1[i], x2[j], x3[k],i,j,k);
	
/* The acceleration vectors are spherical polar */	
	
	gr_out[k][j][i]=g[0];
	gt_out[k][j][i]=g[1];
	gp_out[k][j][i]=g[2];
    
//    dv_ds_out[k][j][i]=disk_dv_ds_array[k][j][i];
    
    

    M_max=4400.;    
 //   eta_max=pow((M_max/(krad*(1.-alpha))),(1./alpha));
    
    eta_max=7.954346e+07;
    
#if EOS!=ISOTHERMAL
    T   = v[PRS]/v[RHO]*KELVIN*mu;    //Compute initial temperature in Kelvin
#else
	T=g_inputParam[T_ISO];               //isothermal temperature
	
#endif
	
	
//    T   = Vc[PRS]/Vc[RHO]*KELVIN*mu;    //Compute initial temperature in Kelvin
    krad=g_inputParam[KRAD];               
    alpharad=g_inputParam[ALPHARAD]; 
	
    if (krad==999 && alpharad==999)
	{
		for (ii=0;ii<MPOINTS;ii++)
		{
			M_UV_array[ii]=M_UV_fit[ii][k][j][i];
		}
	}
	M_temp1=-1.0;
	M_temp2=-1.0;
	fmax=1e-99;
	for (iangle=0;iangle<NFLUX_ANGLES;iangle++)
	{	 
	    if (dvds_array[iangle][k][j][i]>0.0)
	    {				
	        t_UV=   sigma_e * rho * v_th / dvds_array[iangle][k][j][i];  				
		    if (krad==999 && alpharad==999)
			{
				M_UV=linterp(log10(t_UV),t_fit,M_UV_array,MPOINTS);
			}
			else
			{
		        M_UV=krad*pow(t_UV,(alpharad));	
			}				
	        if (M_UV>M_max)
	            M_UV=M_max;
			if (M_UV>M_temp1)
				M_temp1=M_UV;
	    }
        if (pow(flux_r_UV[iangle][k][j][i],2)+pow(flux_t_UV[iangle][k][j][i],2)>fmax)
		{
			fmax=pow(M_UV*flux_r_UV[iangle][k][j][i],2)+pow(M_UV*flux_t_UV[iangle][k][j][i],2);
			M_temp2=M_UV;
		}
		
		
		
	}
	M_out1[k][j][i]=M_temp1;
	M_out2[k][j][i]=M_temp2;
	
	#endif
	
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





