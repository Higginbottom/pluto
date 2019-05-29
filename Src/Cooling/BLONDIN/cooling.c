/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Take a source step using power-law cooling.


 
  \authors Nick Higginbottom
  \date    July 28, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define ITMAX 100
#define EPS 3.0e-8

double xi,sqsqxi,sqxi,ne,nH,n,tx,E,hc_init,rho,dt_share;
double comp_c_pre,comp_h_pre,line_c_pre,brem_c_pre,xray_h_pre;
double heatcool();
double zfunc();
double zbrent();
double zfunc2();

int flag;

/* ********************************************************************* */
void BlondinCooling (const Data *data, double dt, timeStep *Dts, Grid *grid)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int     i, j, k,iii;
  double tmin,tmax,dT;
  double  cost, dE, E1,E_f;
  double  p, T, p_f, T_f;
  double  lx, r, test;
  double t_u,t_l,mu,t_new;
  double dt_min,***xi_out,***T_out;
  double hc_final;
  double T_test;
  int imin,jmin,reset_flag;
  double tstore1,tstore2,nenh_store;
  double ***dt_out;
  dt_out     = GetUserVar("dt");

  dt_share=dt*UNIT_TIME;  //We need to share the current time step so the zbrent code can use it - must be in real units
  lx=g_inputParam[L_x];  //Xray luminosiy
  tx=g_inputParam[T_x];  //Xray tenperature
  mu=g_inputParam[MU];  //Mean particle mass  
  flag=0;
  
  
/*  mu=MeanMolecularWeight(v); This is how it should be done - but we are ionized and get the wrong answer */
  
//  mu=0.6;
  

  reset_flag=0;

  dE = 1.e-18;  //This is from the original cooling code - unsure if I need it.....
  
  /* DOM_LOOP below does a loop from   
  j(theta) from 2 to 101
  i(r) from 2 to 201
    this is for 100 theta and 200 r points - so ignores ghost zones...  
  */
  
  DOM_LOOP(k,j,i){
	  
	comp_c_pre=data->comp_c_pre[k][j][i];
	comp_h_pre=data->comp_h_pre[k][j][i];
	line_c_pre=data->line_c_pre[k][j][i];
	brem_c_pre=data->brem_c_pre[k][j][i];
	xray_h_pre=data->xray_h_pre[k][j][i];
	  
	  


	r=grid->x[IDIR][i]*UNIT_LENGTH;  //The radius - in real units
	
	
	
	
	


    rho = data->Vc[RHO][k][j][i]*UNIT_DENSITY;  //density of the current cell in physical units
	 
    p   = data->Vc[PRS][k][j][i];  //pressure of the current cell
    T   = data->Vc[PRS][k][j][i]/data->Vc[RHO][k][j][i]*KELVIN*mu;    //Compute initial temperature in Kelvin
    E   = (p*UNIT_PRESSURE)/(g_gamma-1);     //Compute internal energy in physical units
	 
	nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming stellar abundances
	
	xi=lx/nH/r/r;     //ionization parameter
	ne=nH*ne_rat(T,xi);
	ne=1.21*nH;
//	ne=1.21*nH;             //electron number density assuming full ionization
	
	n=rho/(mu*CONST_mp);    //particle density

	 
	T=E*(2.0/3.0)/(n*CONST_kB);
	
	
	  
    if (T < g_minCoolingTemp)
		{
			continue;  //Quit if the temperature is too cold - this may need tweeking
		}	
	
	

	
//E1   = T*n*CONST_kB/(g_gamma-1);
	
	
	
	sqsqxi=pow(xi,0.25);    //we use xi^0.25 in the cooling rates - expensive to recompute every call, so do it now and transmit externally	
	sqxi=sqrt(xi);
	hc_init=heatcool(T);    //Get the initial heating/cooling rate


//   the next few lines bracket the solution temperature
	t_l=T*0.9;
	t_u=T*1.1;
	test=zfunc(t_l)*zfunc(t_u);
 	while (test > 0)
	{
		t_l=t_l*0.9;
		t_u=t_u*1.1;
  		test=zfunc(t_l)*zfunc(t_u);
	}
	


    T_f=zbrent(zfunc,t_l,t_u,1.0);
	flag=0;
	
	hc_final=heatcool(T_f);


	if (hc_final*hc_init<0.)   //We have crossed the equilibrium temperature
	{			
		T_test=zbrent(heatcool,fmin(T_f,T),fmax(T_f,T),1.0);
//		if (i==7 && j==86)
//			printf ("We have crossed equilibrium %e %e %e\n",T,T_f,T_test);
		T_f=T_test;
	}

	
	
/*  ----  Update Energy  ----  */
	T_f = MAX (T_f, g_minCoolingTemp); //if the temperature has dropped too low - use the floor (50K)
	
//	heatcool2(data,T,i,j,k);

	E_f=T_f/(2.0/3.0)*(n*CONST_kB); //convert back to energy
	
	p_f=E_f*(g_gamma-1)/UNIT_PRESSURE; 
	
//    p_f = T_f*VV[RHO][k][j][i]/(KELVIN*mu);  //Compute the new pressure from the new temperature - code units
	 	
	
	/*I need to understand this a bit more clearly - p_f/p will give the ratio of the new pressure over the initial pressure
	so 1 minus that will give zero if there is no change - the +1e-18 is there to stop a divide by zero error in the next step
	if that happens. */
	
	
    dE = (fabs(1.0 - E_f/E)) + 1.e-18;  //The fractional change in pressure (or energy since they are proportional)

	
    data->Vc[PRS][k][j][i] = p_f;  //Set the pressure in the cell to the new value

	/* This is a bit obscure - it is from the original code, and appears to set the cooling timescale
	to a value that means the change in energy will be less than the max cooling rate. It needs a bit 
	more understanding - to see if it is what I want. g_maxCoolingRate is the "maximum fractional variation
	due to cooling from one step to the next" and is set to 0.1 in globals.h 
	
	max cooling rate / dE will be equal to 1 if our change in energy is at the max, so dt will equal the largest
	acceptable time step. If this is less than the current time step then we will be reducing it next time. 
	*/
	if (Dts->dt_cool> dt*g_maxCoolingRate/dE)
	{
//		printf ("%i %i %e %e\n",i,j,Dts->dt_cool,dt*g_maxCoolingRate/dE);
		
    Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE); 
	imin=i;
	jmin=j;
	tstore1=T;
	tstore2=T_f;
	nenh_store=ne*nH*ne_rat(T,xi);
	reset_flag=1;
	}
	dt_out[k][j][i]=dt*g_maxCoolingRate/dE;
	
//	if (i==7 && j==86)
//	{
//	printf ("cell1 %i %i xi %e T_i %e T_f %e dt %e dE/dt %e E %e Ef %e P %e Pf %e time %e test %e\n",i,j,xi,T,T_f,dt*g_maxCoolingRate/dE,(E_f-E)/dt,E,E_f,p,p_f,g_time,hc_init*hc_final);
	//} 
//	if (i==7 && j==87)
//	{
//	printf ("cell2 %i %i xi %e T_i %e T_f %e dt %e dE/dt %e E %e Ef %e P %e Pf %e time %e test %e\n",i,j,xi,T,T_f,dt*g_maxCoolingRate/dE,(E_f-E)/dt,E,E_f,p,p_f,g_time,hc_init*hc_final);
	//} 
/*
	
	if (g_time>0.05)

	{
//		exit(0);
		
	flag=1;
	
	tmin=20000.;
	tmax=25000.;
	
	dT=(tmax-tmin)/1000.;
	for (iii=0;iii<1001;iii++)
	{
		heatcool(tmin+iii*dT);
		zfunc(tmin+iii*dT);
	}
	
	exit(0);
	
	
	
	flag=0;
	}
	}*/
  }
//   if (reset_flag==1)
//    {
//  printf ("min_dt i=%i j=%i T_i %e T_f %e nenh %e\n",imin,jmin,tstore1,tstore2,nenh_store);
  //}
//printf ("%e %e %e %e\n",data->Vc[RHO][0][98][3]*UNIT_DENSITY,data->Vc[RHO][0][99][3]*UNIT_DENSITY,data->Vc[RHO][0][100][3]*UNIT_DENSITY,data->Vc[RHO][0][102][3]*UNIT_DENSITY);
}




double heatcool(double T)
{
	double lambda,st,h_comp,c_comp,h_xray,l_line,l_brem;
	ne=nH*ne_rat(T,xi);
	
	st=sqrt(T);
	h_comp=comp_h_pre*8.9e-36*xi*tx*(ne*nH);
	c_comp=comp_c_pre*8.9e-36*xi*(4.0*T)*ne*nH;
	h_xray=xray_h_pre*1.5e-21*(sqsqxi/st)*(1-(T/tx))*nH*nH;
//	l_line=line_c_pre*(1.7e-18*exp(-1.3e5/T)/(xi*st)+1e-24)*ne*nH;	
	
	l_line=(line_c_pre*((1e-16*exp(-1.3e5/T)/sqxi/T)+fmin(fmin(1e-24,5e-27*st),1.5e-17/T)))*ne*nH;
	
	
	l_brem=brem_c_pre*3.3e-27*st*ne*nH;	
	lambda=h_comp+h_xray-l_brem-l_line-c_comp;
	if (flag==1)
	{
		printf ("T %e xi %e h_comp %e c_comp %e h_xray %e l_line %e l_brem %e lambda %e E %e dt %e\n",T,xi,h_comp,c_comp,h_xray,l_line,l_brem,lambda,T*n*CONST_kB/(2.0/3.0),dt_share);
	}

	return (lambda);
	
	
}


//A little copy for diagnostic purposes

double heatcool2(double xi,double T,int i, int j,int k, double ne, double nh)
{
	double lambda,st,h_comp,c_comp,h_xray,l_line,l_brem;
	double ***comp_h, ***comp_c, ***line_c, ***brem_c, ***xray_h ;
	


	comp_h  = GetUserVar("ch");
	comp_c  = GetUserVar("cc");
	xray_h  = GetUserVar("xh");
	line_c  = GetUserVar("lc");
	brem_c  = GetUserVar("bc");
	
	
		
	st=sqrt(T);
	comp_h[k][j][i]=h_comp=8.9e-36*xi*tx;
	comp_c[k][j][i]=c_comp=8.9e-36*xi*(4.0*T);
	xray_h[k][j][i]=h_xray=1.5e-21*(pow(xi,0.25)/st)*(1-(T/tx));
//	line_c[k][j][i]=l_line=(1.7e-18*exp(-1.3e5/T)/(xi*st)+1e-24);
	line_c[k][j][i]=l_line=(1e-16*exp(-1.3e5/T)/sqrt(xi)/T)+fmin(fmin(1e-24,5e-27*st),1.5e-17/T);
		
	brem_c[k][j][i]=l_brem=3.3e-27*st;	
	lambda=(h_comp+h_xray-l_brem-l_line-c_comp);
//	lambda=h_comp-c_comp-l_brem-l_line;
//	printf ("BLAH5 %e %e %e %e %e %e %e %e %e %e %e\n",g_time,xi,T,E,h_comp,c_comp,l_brem,l_line,h_xray,lambda,dt_share);	
	return (lambda);
	
	
}


/*This function calculates 

the internal energy at temperature (temp),
minus the current energy
minus the *average* heating/cooling rate over the time step dt.

this will equal zero if we have the 'correct' final temperature. 

It uses several external variables to allow zbrent to be used to find the root*/


double zfunc(double temp) 
{	
	double ans;
	ans=(temp*n*CONST_kB/(2.0/3.0))-E-dt_share*(hc_init+heatcool(temp))/2.0;
	if (flag==1)
		printf ("t %e ANS= %e hc_init %e heatcool %e E1 %e E2 %e mean_hc %e\n",temp,ans,hc_init,heatcool(temp),E,(temp*n*CONST_kB/(2.0/3.0)),dt_share*(hc_init+heatcool(temp))/2.0);
	return (ans);
}





//A simple copy of zbrent from python....

double
zbrent (func, x1, x2, tol)
     double x1, x2, tol;
     double (*func) (double);   /* ANSI: double (*func)(double); */
{
  int iter;
  double a = x1, b = x2, c, d, e, min1, min2;
  double fa = (*func) (a), fb = (*func) (b), fc, p, q, r, s, tol1, xm;

  c = d = e = 0;                // to avoid -03 warning


  if (fb * fa > 0.0)
  {
    printf ("ZBRENT: Min %e & Max %e must bracket zero, but got %e & %e\n", x1, x2, fa, fb);
  }
  fc = fb;
  for (iter = 1; iter <= ITMAX; iter++)
  {
    if (fb * fc > 0.0)
    {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs (fc) < fabs (fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * EPS * fabs (b) + 0.5 * tol;
    xm = 0.5 * (c - b);
    if (fabs (xm) <= tol1 || fb == 0.0)
      return b;
    if (fabs (e) >= tol1 && fabs (fa) > fabs (fb))
    {
      s = fb / fa;
      if (a == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = fabs (p);
      min1 = 3.0 * xm * q - fabs (tol1 * q);
      min2 = fabs (e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2))
      {
        e = d;
        d = p / q;
      }
      else
      {
        d = xm;
        e = d;
      }
    }
    else
    {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs (d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs (tol1) : -fabs (tol1));
    fb = (*func) (b);
  }
  printf ("Maximum number of iterations exceeded in ZBRENT\n");
  return b;
}

#undef ITMAX
#undef EPS

