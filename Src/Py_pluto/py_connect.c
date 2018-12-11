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
	
	
	printf ("DOING IT\n");
	if (flag==0) //Initialising first time round
    DOM_LOOP(k,j,i){
		d->comp_h_pre[k][j][i]=1.0;
		d->comp_c_pre[k][j][i]=1.0;
		d->xray_h_pre[k][j][i]=1.0;
		d->line_c_pre[k][j][i]=1.0;
		d->brem_c_pre[k][j][i]=1.0;
	}
	else
	{
		if ((fptr = fopen ("prefactors.dat", "r")) == NULL)
		{
			printf ("NO prefactor file\n");
			exit(0);
		}
	    DOM_LOOP(k,j,i){ //Initialise
			d->comp_h_pre[k][j][i]=-1.0;
			d->comp_c_pre[k][j][i]=-1.0;
			d->xray_h_pre[k][j][i]=-1.0;
			d->line_c_pre[k][j][i]=-1.0;
			d->brem_c_pre[k][j][i]=-1.0;
		}
		fgets (aline, LINELENGTH, fptr);
		icount=0;
		while (fgets (aline, LINELENGTH, fptr) != NULL)	
		{		
			if ((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le", &ii, &rcen, &jj, &thetacen, &dens, &comp_h_pre, &comp_c_pre, &xray_h_pre, &brem_c_pre, &line_c_pre)) == 10)
			{
				DOM_LOOP(k,j,i)
				{
					if (fabs(((rcen/UNIT_LENGTH)-grid->x[IDIR][i])/(rcen/UNIT_LENGTH))<1e-6 && fabs((thetacen-grid->x[JDIR][j])/thetacen)<1e-6)
					{
						if ((d->Vc[RHO][k][j][i]*UNIT_DENSITY/dens)-1.>1e-6)
						{
							printf ("Density mismatch\n");
						}
						else
						{
							icount++;
							d->comp_h_pre[k][j][i]=comp_h_pre;
							d->comp_c_pre[k][j][i]=comp_c_pre;
							d->xray_h_pre[k][j][i]=xray_h_pre;
							d->line_c_pre[k][j][i]=line_c_pre;
							d->brem_c_pre[k][j][i]=brem_c_pre;
						}
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
				printf ("PROBLEM with i=%i j=%i\n",i,j);
			}
		}
		printf ("Read in %i prefectors\n",icount);
	}
}



