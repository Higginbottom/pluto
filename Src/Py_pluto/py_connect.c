/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Read in a py_heatcool file and compare with known heating and cooling rates

n constant
  
 
  \authors Nick H
  \date    Dec 3 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define LINELENGTH 400



/* ********************************************************************* */
void read_py_heatcool (Data *d, Grid *grid,int flag)
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
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	int ii,jj,nwords;
	double dens,comp_h_pre,comp_c_pre,xray_h_pre,brem_c_pre,line_c_pre;
	double rcen,thetacen;
	int icount;
	#if (BODY_FORCE & VECTOR)
	double gx_pre,gy_pre,gz_pre;
	#endif
	
	

	if (flag==0) //Initialising first time round
    DOM_LOOP(k,j,i){
		d->comp_h_pre[k][j][i]=1.0;
		d->comp_c_pre[k][j][i]=1.0;
		d->xray_h_pre[k][j][i]=1.0;
		d->line_c_pre[k][j][i]=1.0;
		d->brem_c_pre[k][j][i]=1.0;
		g_rad_force_pre[0][k][j][i]=1.0;
		g_rad_force_pre[1][k][j][i]=1.0;
		g_rad_force_pre[2][k][j][i]=1.0;		
	}
	else
	{
		fptr = fopen ("prefactors.dat", "r");
		if (fptr==NULL)
		{
		    DOM_LOOP(k,j,i){
				d->comp_h_pre[k][j][i]=1.0;
				d->comp_c_pre[k][j][i]=1.0;
				d->xray_h_pre[k][j][i]=1.0;
				d->line_c_pre[k][j][i]=1.0;
				d->brem_c_pre[k][j][i]=1.0;
				g_rad_force_pre[0][k][j][i]=1.0;
				g_rad_force_pre[1][k][j][i]=1.0;
				g_rad_force_pre[2][k][j][i]=1.0;		
			}
			printf ("NO prefactor file\n");
		}
		else
		{
	    DOM_LOOP(k,j,i){ //Initialise
			d->comp_h_pre[k][j][i]=-1.0;
			d->comp_c_pre[k][j][i]=-1.0;
			d->xray_h_pre[k][j][i]=-1.0;
			d->line_c_pre[k][j][i]=-1.0;
			d->brem_c_pre[k][j][i]=-1.0;
			g_rad_force_pre[0][k][j][i]=1.0;
			g_rad_force_pre[1][k][j][i]=1.0;
			g_rad_force_pre[2][k][j][i]=1.0;
		}
		fgets (aline, LINELENGTH, fptr);
		icount=0;
		while (fgets (aline, LINELENGTH, fptr) != NULL)	
		{		
			if 
				#if (BODY_FORCE & VECTOR)
				((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le %le %le %le", &ii, &rcen, &jj, &thetacen, &dens,
				&comp_h_pre, &comp_c_pre, &xray_h_pre, &brem_c_pre, &line_c_pre,
				&gx_pre, &gy_pre, &gz_pre)) == 13)
				#else
				((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le", &ii, &rcen, &jj, &thetacen, &dens,
				&comp_h_pre, &comp_c_pre, &xray_h_pre, &brem_c_pre, &line_c_pre)) == 10)
				#endif
			{
				DOM_LOOP(k,j,i)
				{
					if (fabs(((rcen/UNIT_LENGTH)-grid->x[IDIR][i])/(rcen/UNIT_LENGTH))<1e-5 && fabs((thetacen-grid->x[JDIR][j])/thetacen)<1e-5)
					{
						if ((d->Vc[RHO][k][j][i]*UNIT_DENSITY/dens)-1.>1e-6)
						{
							printf ("Density mismatch in cell i=%i j=%i old=%e new=%e\n",i,j,d->Vc[RHO][k][j][i]*UNIT_DENSITY,dens);
						}
						icount++;
						d->comp_h_pre[k][j][i]=comp_h_pre;
						d->comp_c_pre[k][j][i]=comp_c_pre;
						d->xray_h_pre[k][j][i]=xray_h_pre;
						d->line_c_pre[k][j][i]=line_c_pre;
						d->brem_c_pre[k][j][i]=brem_c_pre;
						#if (BODY_FORCE & VECTOR)
							g_rad_force_pre[0][k][j][i]=gx_pre;
							g_rad_force_pre[1][k][j][i]=gy_pre;
							g_rad_force_pre[2][k][j][i]=gz_pre;
						#endif
					}
				}
			}
			else
			{
				printf ("Prefactor file incorrectly formatted nwords=%i %s\n",nwords,aline);
				exit(0);
			}
		}
	    DOM_LOOP(k,j,i){ //Test
			if (d->comp_h_pre[k][j][i]<0.0)
			{
				printf ("comp_h_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->comp_h_pre[k][j][i]);
			}
			if (d->comp_c_pre[k][j][i]<0.0)
			{
				printf ("comp_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->comp_c_pre[k][j][i]);
			}
			if (d->xray_h_pre[k][j][i]<0.0)
			{
				printf ("xray_h_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->xray_h_pre[k][j][i]);
			}
			if (d->line_c_pre[k][j][i]<0.0)
			{
				printf ("line_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->line_c_pre[k][j][i]);
			}
			if (d->brem_c_pre[k][j][i]<0.0)
			{
				printf ("brem_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->brem_c_pre[k][j][i]);
			}
		}
	
		printf ("Read in %i prefectors\n",icount);
	}
}
}


double ne_rat(double temp,double xi)
{
	 double x1,x2,ne_rat;

	x1=1e-2+pow(10,(-51.59417133+12.27740153*log10(temp)));
	x2=pow(10,(-3.80749689+0.86092628*log10(temp)));
	
	ne_rat=fmin(x1,x2);
	ne_rat=fmin(ne_rat,1.21);
//	printf ("%e %e %e %e\n",temp,x1,x2,ne_rat);
	return (ne_rat);
 }

