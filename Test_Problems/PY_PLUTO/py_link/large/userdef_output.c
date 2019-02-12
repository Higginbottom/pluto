#include "pluto.h"

double heatcool2();


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
  int i, j, k;  
double ***comp_h_pre, ***comp_c_pre, ***line_c_pre, ***brem_c_pre, ***xray_h_pre, ***T_out, ***xi_out, ***ne_out, ***nH_out;
double rho,p,T,xi,ne,nH,n,mu,lx,tx,r;
comp_h_pre  = GetUserVar("ch_pre");
comp_c_pre  = GetUserVar("cc_pre");
xray_h_pre  = GetUserVar("xh_pre");
line_c_pre  = GetUserVar("lc_pre");
brem_c_pre  = GetUserVar("bc_pre");
ne_out	= GetUserVar("ne");
nH_out	= GetUserVar("nh");
T_out     = GetUserVar("T");
xi_out     = GetUserVar("XI");

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
	  
      ne=nH*ne_rat(T,xi);
	  
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





