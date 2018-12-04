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
	
	
	
	
	
	printf ("DOING IT\n");
	if (flag==1) //Initialising first tome round
    DOM_LOOP(k,j,i){
		d->comp_h_pre[k][j][i]=1.0;
		d->comp_c_pre[k][j][i]=1.0;
		d->xray_h_pre[k][j][i]=1.0;
		d->line_c_pre[k][j][i]=1.0;
		d->brem_c_pre[k][j][i]=1.0;
	}
}



