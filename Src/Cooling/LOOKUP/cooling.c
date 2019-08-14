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





double zbrent();
double zfunc2();
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
    int i,j,nx,ny;
	double xi_test;
	double T_test;
	double dt_test;
	

	if (lookup_flag==0)
	{
		read_heatcool("heatcool_lookup.dat");
		printf ("We will read in the lookup table\n");	
		printf ("Finished reading in everything\n");		
		lookup_flag=1;
	}
	else
	{
		printf ("We are in Lookupcooling lookupflag=%i\n",lookup_flag);
		printf ("We have read in the heatcool file\n");

		dt_test=(T_lu[30]-T_lu[22])/100.;
		xi_test=xi_lu[0];

		

		
		
		for (i=0;i<101;i++)
		{
			T_test=T_lu[22]+i*dt_test;
			printf ("TEST %e %e %e\n",T_test,xi_test,gsl_spline2d_eval(spline,xi_test , T_test, xacc, yacc));
		}



		exit(0);

		
			
//	gsl_interp2d_alloc(const gsl_interp2d_type * T, const size_t xsize, const size_t ysize)				
			

	
	}
}



int
	read_heatcool(char* fname)
{
    FILE *fopen (), *fptr;
	char line[1000];
	int n,number,i,j,count;
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
	
	T= gsl_interp2d_bilinear;
	
	xi_lu=malloc(n_xi_lu*sizeof(double));
	T_lu=malloc(n_T_lu*sizeof(double));
	hc_lu=malloc(n_xi_lu*n_T_lu*sizeof(double));
	
	spline = gsl_spline2d_alloc(T, n_xi_lu, n_T_lu);	
	xacc = gsl_interp_accel_alloc();
	yacc = gsl_interp_accel_alloc();	

	
			
	
	

		
	
	count=0;
	for (i=0;i<n_xi_lu;i++)
	{
		for (j=0;j<n_T_lu;j++)
		{	
			count++;
			printf ("%i %i %i %i\n",i,j,count,n_xi_lu*n_T_lu);
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

