/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Take a source step using power-law cooling.

  Integrate the ODE 
  \f[
       dp_{\rm cgs}/dt_{\rm cgs}
      = -(\Gamma - 1) \Lambda(\rho_{\rm cgs}, T_{\rm cgs})
       \qquad{\rm where}\qquad
       \Lambda = \frac{a_{br}}{(\mu m_H)^2} \rho_{\rm cgs}^2 \sqrt{T_{\rm cgs}}
      \,,\qquad
         a_{br} = 2.e-27   c.g.s
  \f]
  which accounts for bremmstrahlung cooling.
  Here the subscript 'cgs' means that the corresponding quantity is given
  in c.g.s units. 
  We denote with \c mu the molecular weight, \c mH the hydrogen mass (in c.g.s).
  
  The previous integration is carried out analytically since the density
  does not change during this step.
  In non-dimensional form:
  \f[
      \frac{dp}{dt} = -{\rm cost} \;\rho  (\rho p)^{\HALF}
  \f]
 
   [notice that since  \c p/rho=T/KELVIN this is equivalent to:
    \c dp/dt=-cost rho^2 (T/KELVIN)^(1/2) ]
 
  The quantity \c cost is determined by transforming the dimensional
  equation into the non-dimensional one.
  If p, rho and t are in code (non-dimensional) units and if \c L_0,
  \c rho_0, and \c V_0 are the unit length, density and velocity,
  then \c cost is found to be:
  \verbatim
                   a_br * (gamma - 1) * L_0 * rho_0
        cost = -------------------------------------------
                sqrt(kB * mu * mH) * kB * mu * mH * V_0^2
  \endverbatim
   where a_{br} = 2.e-27 (in c.g.s), kB is the Boltmann constant
  
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 28, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define ITMAX 100
#define EPS 3.0e-8

double xi,sqsqxi,ne,nH,n,tx,E,hc_init,rho,dt_share;
double heatcool();
double heatcool2();
double zfunc();
double zbrent();
double zfunc2();

/* ********************************************************************* */
void BlondinCooling (Data_Arr VV, double dt, Time_Step *Dts, Grid *grid)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  cost, dE;
  double  p, T, p_f, T_f;
  double  lx, r, test;
  double t_u,t_l,mu,t_new;
  double dt_min;


  dt_share=dt*UNIT_TIME;  //We need to share the current time step so the zbrent code can use it - must be in real units
  lx=g_inputParam[L_x];  //Xray luminosiy
  tx=g_inputParam[T_x];  //Xray tenperature
  
  
/*  mu=MeanMolecularWeight(v); This is how it should be done - but we are ionized and get the wrong answer */
  
  mu=0.6;
  

  dE = 1.e-18;  //This is from the original cooling code - unsure if I need it.....
  
  /* DOM_LOOP below does a loop from   
  j(theta) from 2 to 101
  i(r) from 2 to 201
    this is for 100 theta and 200 r points - so ignores ghost zones...  
  */
//printf ("BLAH\n");
  
  DOM_LOOP(k,j,i){

	r=grid[IDIR].x[i]*UNIT_LENGTH;  //The radius - in real units


    rho = VV[RHO][k][j][i]*UNIT_DENSITY;  //density of the current cell in physical units
	 
    p   = VV[PRS][k][j][i];  //pressure of the current cell
    T   = (VV[PRS][k][j][i]/VV[RHO][k][j][i]*KELVIN*mu);    //Compute initial temperature in Kelvin
    E   = p*UNIT_PRESSURE/(g_gamma-1);     //Compute internal energy in physical units
	 
	  
    if (T < g_minCoolingTemp) continue;  //Quit if the temperature is too cold - this may need tweeking
	
	
	nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming slar abundances
	ne=1.21*nH;             //electron number density assuming full ionization
	n=rho/(mu*CONST_mp);    //particle density
	
	E   = T*n*CONST_kB/(g_gamma-1);
	
	if (i==2 && j==37)
	{	
//	printf ("BLAH %e %e %e\n",T,E,T*n*CONST_kB/(g_gamma-1));
    }
	
	
	xi=lx/nH/r/r;     //ionization parameter
	sqsqxi=pow(xi,0.25);    //we use xi^0.25 in the cooling rates - expensive to recompute every call, so do it noe and transmit externally	
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
//		printf ("test=%e\n",test);
	}
	
	

	
	
//	if (j==101)
//	{
//		printf ("%e %e %e ",T_f,xi,rho);
//		heatcool2(T_f);
//	}
		


//we now search for the solution temperature
	
    T_f=zbrent(zfunc,t_l,t_u,1.0);
	
	
/*  ----  Update Energy  ----  */
	T_f = MAX (T_f, g_minCoolingTemp); //if the temperature has dropped too low - use the floor (50K)
    p_f = T_f*VV[RHO][k][j][i]/(KELVIN*mu);  //Compute the new pressure from the new temperature - code units
	
	
	/*I need to understand this a bit more clearly - p_f/p will give the ratio of the new pressure over the initial pressure
	so 1 minus that will give zero if there is no change - the +1e-18 is there to stop a divide by zero error in the next step
	if that happens. */
	
	
    dE = (fabs(1.0 - p_f/p))/UNIT_TIME + 1.e-18;  //The fractional change in pressure (or energy since they are proportional)
	
    VV[PRS][k][j][i] = p_f;  //Set the pressure in the cell to the new value

	/* This is a bit obscure - it is from the original code, and appears to set the cooling timescale
	to a value that means the change in energy will be less than the max cooling rate. It needs a bit 
	more understanding - to see if it is what I want. g_maxCoolingRate is the "maximum fractional variation
	due to cooling from one step to the next" and is set to 0.1 in globals.h 
	
	max cooling rate / dE will be equal to 1 if our change in energy is at the max, so dt will equal the largest
	acceptable time step. If this is less than the current time step then we will be reducing it next time. 
	*/
//	 printf ("T=%e t_f=%e\n",T,T_f);
//	if (i==2 && j==37)
//	{
//			printf ("T=%e T_f=%e (%e) density=%e dt=%e\n",T,T_f,T_f-T,rho,dt);
//			heatcool2(T_f);
//			if (T_f>2e6) exit(0);
//		}
//			heatcool2(T_f);
//			exit(0);

    Dts->dt_cool = MIN(Dts->dt_cool, dt*g_maxCoolingRate/dE); 
	Dts->dt_cool = 1e6;
//    printf ("cooling dt=%e unit_time=%e\n",dt_min,UNIT_TIME);
//    exit(0);	
  }
//  exit(0);
}





double heatcool(double T)
{
	double lambda,st,h_comp,c_comp,h_xray,l_line,l_brem;
	
	st=sqrt(T);  //We use this a couple of times - so compute it at the start
	h_comp=8.9e-36*xi*tx*ne*nH;  //Compton heating
	c_comp=8.9e-36*xi*(4.0*T)*ne*nH;  //Compton cooling
	h_xray=1.5e-21*(sqsqxi/st)*(1-(T/tx))*ne*nH;  //PI heating
	l_line=1.7e-18*(exp(-1.3e5/T)/(xi*st)+1e-24)*ne*nH;  //line+recomb cooling
	l_brem=3.3e-27*st*ne*nH;  //Brem cooling
	
	
	lambda=h_comp+h_xray-l_brem-l_line-c_comp;  //Add em all up

	return (lambda);
	
	
}


//A little copy for diagnostic purposes

double heatcool2(double T)
{
	double lambda,st,h_comp,c_comp,h_xray,l_line,l_brem;
	
	st=sqrt(T);
	h_comp=8.9e-36*xi*tx*(ne*nH);
	c_comp=8.9e-36*xi*(4.0*T)*ne*nH;
	h_xray=1.5e-21*(sqsqxi/st)*(1-(T/tx))*nH*nH;
	l_line=1.7e-18*(exp(-1.3e5/T)/(xi*st)+1e-24)*ne*nH;
	l_brem=3.3e-27*st*ne*nH;	
	lambda=h_comp+h_xray-l_brem-l_line-c_comp;
//	lambda=h_comp-c_comp-l_brem-l_line;
	printf ("%e %e %e %e %e\n",h_comp/(ne*nH),c_comp/(ne*nH),l_brem/(ne*nH),l_line/(ne*nH),h_xray/(nH*nH));	
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
	ans=(temp*n*CONST_kB/(g_gamma-1))-E-dt_share*(hc_init+heatcool(temp))/2.0;
	return (ans);
}


double zfunc2(double temp) 
{
	double ans;
	ans=(temp*n*CONST_kB/(g_gamma-1))-E-dt_share*(hc_init+heatcool2(temp))/2.0;
	printf ("Final E=%e Current E=%e initial_hc=%e new_hc=%e dt=%e\n",temp*n*CONST_kB/(g_gamma-1),E,hc_init,heatcool2(temp),dt_share);
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

