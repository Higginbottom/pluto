#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


int lookup_flag;
//double xi_lookup[60];
//double T_lookup[1000];
int n_xi_lu,n_T_lu;
//double lookup[60][1000];

double *xi_lu;
double *T_lu;
double *hc_lu;




gsl_spline2d *spline;
gsl_interp_accel *xacc;
gsl_interp_accel *yacc;

const gsl_interp2d_type *interpolator;