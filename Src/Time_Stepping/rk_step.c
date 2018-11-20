/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Advance equations with Runge Kutta time integrators.

  Main driver for RK split/unsplit integrations and finite difference
  methods (RK3).
  Time stepping include Euler, RK2 and RK3.

  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Nov 21, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* Weight factor for 2nd stage of RK integrators */

#if TIME_STEPPING == RK2
 #define w0  0.5
 #define wc  0.5
#elif TIME_STEPPING == RK3
 #define w0 0.75
 #define wc 0.25
#endif

/* ********************************************************************* */
int AdvanceStep (Data *d, Riemann_Solver *Riemann, 
                 timeStep *Dts, Grid *grid)
/*!
 * Advance the equations by a single time step using unsplit 
 * integrators based on the method of lines.
 *
 * \param [in,out]      d  pointer to Data structure
 * \param [in]    Riemann  pointer to a Riemann solver function
 * \param [in,out]    Dts  pointer to time step structure
 * \param [in]       grid  pointer to array of Grid structures
 *    
 *********************************************************************** */
{
  int  i, j, k, nv;
  static double  one_third = 1.0/3.0;
  static Data_Arr U0, Bs0;
  RBox   box;
if (g_stepNumber>NSHTIME1)   printf ("ADV\n");
 if (g_stepNumber>NSHTIME1)  printf ("IN ADV        %e density  %e pressure %20.15e temp %e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],KELVIN*0.6*d->Vc[PRS][0][NSHj][NSHi]/d->Vc[RHO][0][NSHj][NSHi]);

	static Data_Arr Uhalf;
  
  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

#if (defined PARTICLES) && DIMENSIONAL_SPLITTING == YES
  print ("! AdvanceStep(): particles require DIMENSIONAL_SPLITTING == NO\n");
  QUIT_PLUTO(1);
#endif

/* --------------------------------------------------------
   0. Allocate memory 
   -------------------------------------------------------- */

  if (U0 == NULL){
    U0 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
Uhalf = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);

    #ifdef STAGGERED_MHD
     Bs0 = ARRAY_4D(DIMENSIONS, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
  }

/* --------------------------------------------------------
   1. Predictor step (EULER, RK2, RK3, MH (split), 
                      ChTr (split))

      After baoundaries have been set we flag zones lying 
      in a shock. 
      This is useful for shock flattening or 
      entropy/energy selective update.

      Note: when using FARGO, boundary condition must be 
      set on the *total* velocity while the update step is 
      performed on the *residual* velocity.
      The addition and subtraction operations are 
      automatically performed in Boundary() function.
   -------------------------------------------------------- */

/* -- 1a. Set boundary conditions -- */
if (g_stepNumber>NSHTIME1)   printf ("ADV B4 BC     %e density  %e pressure %20.15e temp %e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],KELVIN*0.6*d->Vc[PRS][0][NSHj][NSHi]/d->Vc[RHO][0][NSHj][NSHi]);

  g_intStage = 1;  
  Boundary (d, ALL_DIR, grid);
  #if (SHOCK_FLATTENING == MULTID) || (ENTROPY_SWITCH) 
  FlagShock (d, grid);
  #endif


if (g_stepNumber>NSHTIME1)   printf ("ADV AF BC     %e density  %e pressure %20.15e temp %e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],KELVIN*0.6*d->Vc[PRS][0][NSHj][NSHi]/d->Vc[RHO][0][NSHj][NSHi]);

if (g_stepNumber>NSHTIME1) printf ("Start %e  density[i-1][j-1]  %e density[i][j-1]  %e density[i+1][j-1]  %e\n",g_time,d->Vc[RHO][0][NSHj-1][NSHi-1],d->Vc[RHO][0][NSHj-1][NSHi],d->Vc[RHO][0][NSHj-1][NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  density[i-1][j]    %e density[i-1][j]  %e density[i+1][j]    %e\n",g_time,d->Vc[RHO][0][NSHj]  [NSHi-1],d->Vc[RHO][0][NSHj]  [NSHi],d->Vc[RHO][0][NSHj]  [NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  density[i-1][j+1]  %e density[i][j+1]  %e density[i+1][j+1]  %e\n",g_time,d->Vc[RHO][0][NSHj+1][NSHi-1],d->Vc[RHO][0][NSHj+1][NSHi],d->Vc[RHO][0][NSHj+1][NSHi+1]);

if (g_stepNumber>NSHTIME1) printf ("Start %e  pressure[i-1][j-1] %e pressure[i][j-1] %e pressure[i+1][j-1] %e\n",g_time,d->Vc[PRS][0][NSHj-1][NSHi-1],d->Vc[PRS][0][NSHj-1][NSHi],d->Vc[PRS][0][NSHj-1][NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  pressure[i-1][j]   %e pressure[i-1][j] %e pressure[i+1][j]   %e\n",g_time,d->Vc[PRS][0][NSHj]  [NSHi-1],d->Vc[PRS][0][NSHj]  [NSHi],d->Vc[PRS][0][NSHj]  [NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  pressure[i-1][j+1] %e pressure[i][j+1] %e pressure[i+1][j+1] %e\n",g_time,d->Vc[PRS][0][NSHj+1][NSHi-1],d->Vc[PRS][0][NSHj+1][NSHi],d->Vc[PRS][0][NSHj+1][NSHi+1]);

if (g_stepNumber>NSHTIME1) printf ("Start %e  V1[i-1][j-1] %e V1[i][j-1] %e V1[i+1][j-1] %e\n",g_time,d->Vc[VX1][0][NSHj-1][NSHi-1],d->Vc[VX1][0][NSHj-1][NSHi],d->Vc[VX1][0][NSHj-1][NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V1[i-1][j]   %e V1[i-1][j] %e V1[i+1][j]   %e\n",g_time,d->Vc[VX1][0][NSHj]  [NSHi-1],d->Vc[VX1][0][NSHj]  [NSHi],d->Vc[VX1][0][NSHj]  [NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V1[i-1][j+1] %e V1[i][j+1] %e V1[i+1][j+1] %e\n",g_time,d->Vc[VX1][0][NSHj+1][NSHi-1],d->Vc[VX1][0][NSHj+1][NSHi],d->Vc[VX1][0][NSHj+1][NSHi+1]);


if (g_stepNumber>NSHTIME1) printf ("Start %e  V2[i-1][j-1] %e V2[i][j-1] %e V2[i+1][j-1] %e\n",g_time,d->Vc[VX2][0][NSHj-1][NSHi-1],d->Vc[VX2][0][NSHj-1][NSHi],d->Vc[VX2][0][NSHj-1][NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V2[i-1][j]   %e V2[i-1][j] %e V2[i+1][j]   %e\n",g_time,d->Vc[VX2][0][NSHj]  [NSHi-1],d->Vc[VX2][0][NSHj]  [NSHi],d->Vc[VX2][0][NSHj]  [NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V2[i-1][j+1] %e V2[i][j+1] %e V2[i+1][j+1] %e\n",g_time,d->Vc[VX2][0][NSHj+1][NSHi-1],d->Vc[VX2][0][NSHj+1][NSHi],d->Vc[VX2][0][NSHj+1][NSHi+1]);

if (g_stepNumber>NSHTIME1) printf ("Start %e  V3[i-1][j-1] %e V3[i][j-1] %e V3[i+1][j-1] %e\n",g_time,d->Vc[VX3][0][NSHj-1][NSHi-1],d->Vc[VX3][0][NSHj-1][NSHi],d->Vc[VX3][0][NSHj-1][NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V3[i-1][j]   %e V3[i-1][j] %e V3[i+1][j]   %e\n",g_time,d->Vc[VX3][0][NSHj]  [NSHi-1],d->Vc[VX3][0][NSHj]  [NSHi],d->Vc[VX3][0][NSHj]  [NSHi+1]);
if (g_stepNumber>NSHTIME1) printf ("Start %e  V3[i-1][j+1] %e V3[i][j+1] %e V3[i+1][j+1] %e\n",g_time,d->Vc[VX3][0][NSHj+1][NSHi-1],d->Vc[VX3][0][NSHj+1][NSHi],d->Vc[VX3][0][NSHj+1][NSHi+1]);



/* -- 1b. Convert primitive to conservative, save initial stage  -- */


 if (g_stepNumber>NSHTIME1)  printf ("ADV B4 PTC     %e density  %e pressure %20.15e temp %e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],KELVIN*0.6*d->Vc[PRS][0][NSHj][NSHi]/d->Vc[RHO][0][NSHj][NSHi]);


   PrimToCons3D(d->Vc, d->Uc, &box);
if (g_stepNumber>NSHTIME1)    printf ("ADV AF PTC    %e density  %20.15e pressure %20.15e current E %20.15e KE=%e INTE=%e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],d->Uc[0][NSHj][NSHi][ENG],0.5*d->Vc[RHO][0][NSHj][NSHi]*(d->Vc[VX1][0][NSHj][NSHi]*d->Vc[VX1][0][NSHj][NSHi]+d->Vc[VX2][0][NSHj][NSHi]*d->Vc[VX2][0][NSHj][NSHi]+d->Vc[VX3][0][NSHj][NSHi]*d->Vc[VX3][0][NSHj][NSHi]),d->Vc[PRS][0][NSHj][NSHi]/(g_gamma - 1.0));
   
  KDOM_LOOP(k) JDOM_LOOP(j){
    memcpy ((void *)U0[k][j][IBEG], d->Uc[k][j][IBEG], NX1*NVAR*sizeof(double));
  }



/* -- 1d. Advance conservative variables array -- */
if (g_stepNumber>NSHTIME1)   printf ("ADV B4 UDS    %e density  %20.15e pressure %20.15e current E %e KE=%e INTE=%e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],d->Uc[0][NSHj][NSHi][ENG],0.5*d->Vc[RHO][0][NSHj][NSHi]*(d->Vc[VX1][0][NSHj][NSHi]*d->Vc[VX1][0][NSHj][NSHi]+d->Vc[VX2][0][NSHj][NSHi]*d->Vc[VX2][0][NSHj][NSHi]+d->Vc[VX3][0][NSHj][NSHi]*d->Vc[VX3][0][NSHj][NSHi]),d->Vc[PRS][0][NSHj][NSHi]/(g_gamma - 1.0));
  UpdateStage(d, d->Uc, NULL, Riemann, g_dt, Dts, grid);
if (g_stepNumber>NSHTIME1)   printf ("ADV AF UDS    %e density  %20.15e pressure %20.15e current E %e KE=%e INTE=%e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],d->Uc[0][NSHj][NSHi][ENG],0.5*d->Vc[RHO][0][NSHj][NSHi]*(d->Vc[VX1][0][NSHj][NSHi]*d->Vc[VX1][0][NSHj][NSHi]+d->Vc[VX2][0][NSHj][NSHi]*d->Vc[VX2][0][NSHj][NSHi]+d->Vc[VX3][0][NSHj][NSHi]*d->Vc[VX3][0][NSHj][NSHi]),d->Vc[PRS][0][NSHj][NSHi]/(g_gamma - 1.0));
  



 
/* -- 1f. Convert to primitive vars -- */
if (g_stepNumber>NSHTIME1)   printf ("ADV B4 CTP    %e density  %e pressure %20.15e current E %e KE=%e INTE=%e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],d->Uc[0][NSHj][NSHi][ENG],0.5*d->Vc[RHO][0][NSHj][NSHi]*(d->Vc[VX1][0][NSHj][NSHi]*d->Vc[VX1][0][NSHj][NSHi]+d->Vc[VX2][0][NSHj][NSHi]*d->Vc[VX2][0][NSHj][NSHi]+d->Vc[VX3][0][NSHj][NSHi]*d->Vc[VX3][0][NSHj][NSHi]),d->Vc[PRS][0][NSHj][NSHi]/(g_gamma - 1.0));
  ConsToPrim3D (d->Uc, d->Vc, d->flag, &box);
if (g_stepNumber>NSHTIME1)   printf ("ADV AF CTP    %e density  %e pressure %20.15e current E %e KE=%e INTE=%e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],d->Uc[0][NSHj][NSHi][ENG],0.5*d->Vc[RHO][0][NSHj][NSHi]*(d->Vc[VX1][0][NSHj][NSHi]*d->Vc[VX1][0][NSHj][NSHi]+d->Vc[VX2][0][NSHj][NSHi]*d->Vc[VX2][0][NSHj][NSHi]+d->Vc[VX3][0][NSHj][NSHi]*d->Vc[VX3][0][NSHj][NSHi]),d->Vc[PRS][0][NSHj][NSHi]/(g_gamma - 1.0));

 if (g_stepNumber>NSHTIME1)  printf ("ADV end       %e density  %e pressure %20.15e temp %e\n",g_time,d->Vc[RHO][0][NSHj][NSHi],d->Vc[PRS][0][NSHj][NSHi],KELVIN*0.6*d->Vc[PRS][0][NSHj][NSHi]/d->Vc[RHO][0][NSHj][NSHi]);



  return 0; /* -- step has been achieved, return success -- */
}

