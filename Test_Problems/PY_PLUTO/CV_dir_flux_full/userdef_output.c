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
  int i, j, k, nv;  
#if COOLING != NO 
double ***comp_h_pre, ***comp_c_pre, ***line_c_pre, ***brem_c_pre, ***xray_h_pre;
#endif

double ***T_out, ***xi_out, ***ne_out, ***nH_out,krad,alpha,krad_UV,alpha_UV,t_UV,M_UV;

#if (BODY_FORCE & VECTOR)
double ***gr_out, ***gt_out, ***gp_out, g[3], Vc[NVAR];
double ***gx_pre, ***gy_pre, ***gz_pre;
#endif

double ***dv_ds_out;
double ***t_out,***M_out,M_max,eta_max,tau_max;

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
M_out=GetUserVar("M");





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
    
    
    krad=0.2;
    alpha=0.6;
    M_max=4400.;    
 //   eta_max=pow((M_max/(krad*(1.-alpha))),(1./alpha));
    
    eta_max=7.954346e+07;
    
 
    
//    if (disk_dv_ds_array[k][j][i]>0.0)
//    {
//        t_out[k][j][i]=   sigma_e * rho * v_th / disk_dv_ds_array[k][j][i];  
//        tau_max=t_UV*eta_max;
//        M_UV=krad_UV*pow(t_UV,(-1.0*alpha_UV))*(pow((1.+tau_max),(1.-alpha_UV))-1.)/(pow(tau_max,(1.-alpha_UV)));
//        M_out[k][j][i]=krad_UV*pow(t_out[k][j][i],(alpha_UV));
//        if (M_out[k][j][i]>M_max)
//            M_out[k][j][i]=M_max;
//    }
//    else
//    {
//        M_out[k][j][i]=M_max; //This is a bit of a fudge - the only way this should be negtaive is if there is no flux..
//    }
	
	
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





