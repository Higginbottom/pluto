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
double ***comp_h_pre, ***comp_c_pre, ***line_c_pre, ***brem_c_pre, ***xray_h_pre, ***T_out, ***xi_out, ***ne_out, ***nH_out;

#if (BODY_FORCE & VECTOR)
double ***gr_out, ***gt_out, ***gp_out, g[3], Vc[NVAR];
double ***gx_pre, ***gy_pre, ***gz_pre;
#endif

double rho,p,T,xi,ne,nH,n,mu,lx,tx,r;
double *x1   = grid->x[IDIR],  *x2  = grid->x[JDIR],  *x3  = grid->x[KDIR];

comp_h_pre  = GetUserVar("ch_pre");
comp_c_pre  = GetUserVar("cc_pre");
xray_h_pre  = GetUserVar("xh_pre");
line_c_pre  = GetUserVar("lc_pre");
brem_c_pre  = GetUserVar("bc_pre");
ne_out	= GetUserVar("ne");
nH_out	= GetUserVar("nh");
T_out     = GetUserVar("T");
xi_out     = GetUserVar("XI");

#if (BODY_FORCE & VECTOR)
gr_out     = GetUserVar("gr");
gt_out     = GetUserVar("gt");
gp_out     = GetUserVar("gp");
gx_pre     = GetUserVar("gx_pre");
gy_pre     = GetUserVar("gy_pre");
gz_pre     = GetUserVar("gz_pre");
#endif


mu=g_inputParam[MU];
lx=g_inputParam[L_x];  //Xray luminosiy
tx=g_inputParam[T_x];  //Xray tenperature 

  
  
  
  DOM_LOOP(k,j,i){
  	  r=grid->x[IDIR][i]*UNIT_LENGTH;  //The radius - in real units
	  
	  
      rho = d->Vc[RHO][k][j][i]*UNIT_DENSITY;  //density of the current cell in physical units	 
      p   = d->Vc[PRS][k][j][i];  //pressure of the current cell
      T   = T_out[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*mu;    //Compute initial temperature in Kelvin
	  

	  
	  nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming stellar abundances
	  
  	  xi_out[k][j][i]=xi=lx/nH/r/r;     //ionization parameter
#if PY_CONNECT
      ne=nH*ne_rat(T,xi);
#else
	  ne=nH*1.21;
#endif
	  
//	  ne=1.21*nH;             //electron number density assuming full ionization
	  n=rho/(mu*CONST_mp);    //particle density
	  
      ne_out[k][j][i]=ne;
      nH_out[k][j][i]=nH; 
	  
	  
	  
	 
	 heatcool2(xi,T,i,j,k,ne,nH);
	  
	comp_h_pre[k][j][i]=d->comp_h_pre[k][j][i];
	comp_c_pre[k][j][i]=d->comp_c_pre[k][j][i];
	xray_h_pre[k][j][i]=d->xray_h_pre[k][j][i];
	line_c_pre[k][j][i]=d->line_c_pre[k][j][i];
	brem_c_pre[k][j][i]=d->brem_c_pre[k][j][i];
	
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
	
/* The prefactors are cartesian */	
	
	gx_pre[k][j][i]=g_rad_force_pre[0][k][j][i];
	gy_pre[k][j][i]=g_rad_force_pre[1][k][j][i];
	gz_pre[k][j][i]=g_rad_force_pre[2][k][j][i];
	
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





