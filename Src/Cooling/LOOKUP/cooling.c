/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Take a source step using power-law cooling.


 
  \authors Nick Higginbottom
  \date    July 28, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define ITMAX 100
#define EPS 3.0e-8


double xi,ne,nH,dt_share,hc_init,E,n;

double heatcool();
double zbrent();
double zfunc();
int read_heatcool(char*);

int flag;

/* ********************************************************************* */
void LookupCooling (Data_Arr VV,const Data *data, double dt, timeStep *Dts, Grid *grid)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
    int i,j,k;
	double T,T_f,t_u,t_l,T_test,test;
	double lx,mu,rho;
	double r,p,p_f,E_f,hc_final;
	
	
    dt_share=dt*UNIT_TIME;  //We need to share the current time step so the zbrent code can use it - must be in real units
    lx=g_inputParam[L_x];  //Xray luminosiy
    mu=g_inputParam[MU];  //Mean particle mass  

	if (lookup_flag==0) //If this is the first time through, set up the interpolators
	{
		read_heatcool("heatcool_lookup.dat");		
		lookup_flag=1;
	}
	
    DOM_LOOP(k,j,i){
	
		r=grid->x[IDIR][i]*UNIT_LENGTH;  //The radius - in real units
	    rho = VV[RHO][k][j][i]*UNIT_DENSITY; 
	    p   = VV[PRS][k][j][i];  //pressure of the current cell
	    T   = VV[PRS][k][j][i]/VV[RHO][k][j][i]*KELVIN*mu;    //Compute initial temperature in Kelvin
	    E   = (p*UNIT_PRESSURE)/(g_gamma-1);     //Compute current internal energy in physical units	 
		nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming stellar abundances	
		xi=lx/nH/r/r;     //ionization parameter
		ne=1.21*nH;             //electron number density assuming full ionization	
		n=rho/(mu*CONST_mp);    //particle density
	 
		T=E*(2.0/3.0)/(n*CONST_kB);
	
		if (T < T_lu[0]) 
		{
//			printf ("T below lowest lookup T %e %e\n",T,T_lu[0]);
			continue;  //Quit if the temperature is too cold - this may need tweeking		
		}
		else if (T > T_lu[n_T_lu-1]) 
		{
//			printf ("T above highest lookup T %e %e\n",T,T_lu[n_T_lu-1]);
			continue;  //Quit if the temperature is too cold - this may need tweeking		
		}
		else if (xi > xi_lu[n_xi_lu-1]) 
		{
//			printf ("xi above highest lookup xi %e %e\n",xi,xi_lu[n_xi_lu-1]);
			continue;  //Quit if the temperature is too cold - this may need tweeking		
		}
		else if (xi < xi_lu[0]) 
		{
//			printf ("xi below lowest lookup xi %e %e\n",xi,xi_lu[0]);
			continue;  //Quit if the temperature is too cold - this may need tweeking		
		}
		hc_init=heatcool(T);    //Get the initial heating/cooling rate
	
//   the next few lines bracket the solution temperature
	t_l=T*0.9;
	t_u=T*1.1;
	
	test=zfunc(t_l)*zfunc(t_u);

 	while (test > 0 && test==test)
	{
		t_l=t_l*0.9;
		t_u=t_u*1.1;
  		test=zfunc(t_l)*zfunc(t_u);
	}
	if (test!=test)  //Test has returned a NAN - which means T cannot be bracketed
	{
		T_f=T;
	}
	else  //We are fine - look for a solution
	{
    T_f=zbrent(zfunc,t_l,t_u,1.0);
	 hc_final=heatcool(T_f);

	if (hc_final*hc_init<0.)   //We have crossed the equilibrium temperature
	{			
		T_test=zbrent(heatcool,fmin(T_f,T),fmax(T_f,T),1.0); //Find the equilibrium
		T_f=T_test;
	}

	}
	
/*  ----  Update Energy  ----  */
	T_f = MAX (T_f, g_minCoolingTemp); //if the temperature has dropped too low - use the floor (50K)
	
//	heatcool2(data,T,i,j,k);

	E_f=T_f/(2.0/3.0)*(n*CONST_kB); //convert back to energy
	
	p_f=E_f*(g_gamma-1)/UNIT_PRESSURE; //and pressure
	
    VV[PRS][k][j][i] = p_f;  //Set the pressure in the cell to the new value
 }
}


int
	read_heatcool(char* fname)
{
    FILE *fopen (), *fptr;
	char line[1000];
	int n,number,i,j;
	double x,y,z;
	char label[1000];
	printf ("We will read in the lookup table\n");
	
    if ((fptr = fopen (fname, "r")) == NULL)
    {
      printf ("No heatcool_lookup.dat file\n");
      exit(0);
    }
	fgets (line, 1000, fptr);
    n = sscanf (line, "%s %d", label, &number);
	if (strncmp (label, "n_xi", 4) == 0) 
	{
		printf ("We have %i xi points\n",number);
		n_xi_lu=number;
	}
	fgets (line, 1000, fptr);
    n = sscanf (line, "%s %d", label, &number);
	if (strncmp (label, "n_t", 3) == 0) 
	{
		printf ("We have %i t points\n",number);
		n_T_lu=number;
	}
	
	interpolator= gsl_interp2d_bilinear;
	
	xi_lu=malloc(n_xi_lu*sizeof(double));
	T_lu=malloc(n_T_lu*sizeof(double));
	hc_lu=malloc(n_xi_lu*n_T_lu*sizeof(double));
	
	spline = gsl_spline2d_alloc(interpolator, n_xi_lu, n_T_lu);	
	xacc = gsl_interp_accel_alloc();
	yacc = gsl_interp_accel_alloc();		
	
	for (i=0;i<n_xi_lu;i++)
	{
		for (j=0;j<n_T_lu;j++)
		{	
			fgets (line, 1000, fptr);
			n = sscanf (line, "%lf %lf %lf",&x,&y,&z);
			xi_lu[i]=x;
			T_lu[j]=y;
			gsl_spline2d_set(spline,hc_lu,i,j,z);
		}
	}
	
	gsl_spline2d_init(spline,  xi_lu, T_lu,hc_lu, n_xi_lu, n_T_lu);

		                
	
	
	
	
	printf ("Finished reading in everything\n");
	
	return(0);
}


double heatcool(double T)
{	
	double lambda,T_temp;
	if (T < T_lu[0]) 
	{
		T=T_lu[0];
	}
	else if (T > T_lu[n_T_lu-1]) 
	{
		T=T_lu[n_T_lu-1];
	}
	lambda=gsl_spline2d_eval(spline,xi,T, xacc, yacc)*nH*ne;
	return (lambda);
}



double heatcool2(double xi,double T,int i, int j,int k, double ne, double nh)
{
	double lambda,st,h_comp,c_comp,h_xray,l_line,l_brem;
	double ***comp_h, ***comp_c, ***line_c, ***brem_c, ***xray_h ;
	
	if (lookup_flag==0) //If this is the first time through, set up the interpolators
	{
		read_heatcool("heatcool_lookup.dat");		
		lookup_flag=1;
	}

	comp_h  = GetUserVar("ch");
	comp_c  = GetUserVar("cc");
	xray_h  = GetUserVar("xh");
	line_c  = GetUserVar("lc");
	brem_c  = GetUserVar("bc");
	
	
	
	comp_c[k][j][i]=0.0;
	xray_h[k][j][i]=0.0;
	line_c[k][j][i]=0.0;
	brem_c[k][j][i]=0.0;
	
	
		
	if (T < T_lu[0]) 
	{
		comp_h[k][j][i]=-1.;
	}
	else if (T > T_lu[n_T_lu-1]) 
	{
		comp_h[k][j][i]=-1.;
	}
	lambda=comp_h[k][j][i]=gsl_spline2d_eval(spline,xi,T, xacc, yacc);

	return (lambda);
	
	
}




double zfunc(double temp) 
{	
	double ans;
	ans=(temp*n*CONST_kB/(2.0/3.0))-E-dt_share*(hc_init+heatcool(temp))/2.0;
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

