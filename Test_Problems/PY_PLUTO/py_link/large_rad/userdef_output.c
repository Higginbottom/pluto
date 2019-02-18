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
  int i, j, k, nv;  
double ***comp_h_pre, ***comp_c_pre, ***line_c_pre, ***brem_c_pre, ***xray_h_pre, ***T_out, ***xi_out, ***ne_out, ***nH_out;

#if (BODY_FORCE & VECTOR)
double ***g1_out, ***g2_out, ***g3_out, g[3], Vc[NVAR];
double ***g1_pre, ***g2_pre, ***g3_pre;
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
g1_out     = GetUserVar("g1");
g2_out     = GetUserVar("g2");
g3_out     = GetUserVar("g3");
g1_pre     = GetUserVar("g1_pre");
g2_pre     = GetUserVar("g2_pre");
g3_pre     = GetUserVar("g3_pre");
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
	
	#if (BODY_FORCE & VECTOR)
	NVAR_LOOP(nv)  Vc[nv]=d->Vc[nv][k][j][i];
	
	BodyForceVector(Vc, g, x1[i], x2[j], x3[k],i,j,k); //This returns r,theta,phi - we need w,y,z
	g1_out[k][j][i]=g[0]*sin(x2[j])+g[1]*cos(x2[j]);
	g2_out[k][j][i]=g[1];		
	g3_out[k][j][i]=g[0]*cos(x2[j])-g[1]*sin(x2[j]);
	
	g1_pre[k][j][i]=g_rad_force_pre[0][k][j][i];
	g2_pre[k][j][i]=g_rad_force_pre[1][k][j][i];
	g3_pre[k][j][i]=g_rad_force_pre[2][k][j][i];
	
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





