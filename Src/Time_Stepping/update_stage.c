/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief   Single stage integration for RK time stepping.

  Advance the equations in conservative form by taking a single stage 
  in the form 
  \f[
      U \quad\Longrightarrow \quad  U + \Delta t R(V)
  \f]
  where \c U is a 3D array of conservative variables, \c V is a 3D array
  of primitive variables, \c R(V) is the right hand side containing 
  flux differences and source terms.
  Note that \c U and \c V may \e not necessarily be the map of 
  each other, i.e., \c U is \e not \c U(V).
  The right hand side can contain contributions from 
   
    - the direction set by the global variable ::g_dir, 
      when DIMENSIONAL_SPLITTING == YES;
    - all directions when DIMENSIONAL_SPLITTING == NO;
   
  When the integrator stage is the first one (predictor), this function 
  also computes the maximum of inverse time steps for hyperbolic and 
  parabolic terms (if the latters are included explicitly).
  
  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
           T. Matsakos

  \date   Sep 07, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void UpdateStage(Data *d, Data_Arr UU, double **aflux,
                 Riemann_Solver *Riemann, double dt, timeStep *Dts, 
                 Grid *grid)
/*!
 * 
 * \param [in,out]  d        pointer to PLUTO Data structure
 * \param [in,out]  UU       data array containing conservative variables
 *                           at the previous time step to be updated
 * \param [out]     aflux    interface fluxes needed for refluxing operations 
 *                           (only with AMR)
 * \param [in]      Riemann  pointer to a Riemann solver function
 * \param [in]      dt       the time step for the current update step
 * \param [in,out]  Dts      pointer to time step structure
 * \param [in]      grid     pointer to Grid structure
 *********************************************************************** */
{
  int  i, j, k;
  int  nv, dir, beg_dir, end_dir;
  int  *ip;

  static Sweep sweep;
  State *stateC = &(sweep.stateC);
  State *stateL = &(sweep.stateL);

  double *inv_dl, dl2;
  static double ***C_dt;
  RBox  sweepBox;

  printf ("In UDS  %e current E %e \n",g_time,UU[0][2][2][ENG]);


#if DIMENSIONAL_SPLITTING == YES
  beg_dir = end_dir = g_dir;
#else
  beg_dir = 0;
  end_dir = DIMENSIONS-1;
#endif

/* --------------------------------------------------------
   0. Allocate memory & reset arrays.
      C_dt is an array used to store the inverse time
      step for the hyperbolic solve.
   -------------------------------------------------------- */

  if (stateC->v == NULL){
    MakeState (&sweep);
    #if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
    C_dt = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #endif
  }

#if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
  if (g_intStage == 1){
    KTOT_LOOP(k) JTOT_LOOP(j){
      memset ((void *)C_dt[k][j],'\0', NX1_TOT*sizeof(double));
    }
  }
#endif




/* --------------------------------------------------------
   2. Update conservative solution array with hyperbolic 
      terms only.
   -------------------------------------------------------- */

  

  for (dir = beg_dir; dir <= end_dir; dir++){

    g_dir = dir;  
	printf ("DIR=%i \n",g_dir);
  /* -- 2b. Set integration box for current update -- */

    RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
    RBoxSetDirections (&sweepBox, g_dir);
    SetVectorIndices (g_dir);


    ResetState(d, &sweep, grid);

    int ntot = grid->np_tot[g_dir];
    int nbeg = *sweepBox.nbeg;
    int nend = *sweepBox.nend;
  
    BOX_TRANSVERSE_LOOP(&sweepBox, k,j,i){
      ip  = sweepBox.n;
      g_i = i;  g_j = j;  g_k = k;
      for ((*ip) = 0; (*ip) < ntot; (*ip)++) {
        NVAR_LOOP(nv) stateC->v[*ip][nv] = d->Vc[nv][k][j][i];
        sweep.flag[*ip] = d->flag[k][j][i];
        #ifdef STAGGERED_MHD
        sweep.bn[*ip] = d->Vs[g_dir][k][j][i];
        #endif
      }


      
      CheckNaN (stateC->v, 0, ntot-1,0);
		if (i==2 || j==2) printf ("Going to States\n");
      States  (&sweep, nbeg - 1, nend + 1, grid);
		if (i==2 || j==2) printf ("Back from States\n");
		
	  if (i==2 || j==2) printf ("Off to riemann i=%i j=%i k=%i beg=%i end=%i\n",i,j,k,nbeg-1,nend);
      Riemann (&sweep, nbeg - 1, nend, Dts->cmax, grid);
	  if (i==2 || j==2) printf ("Back from riemann cmax=%e\n",Dts->cmax);
	  
	  
	  
      if (i==2 || j==2) {
		  printf ("In UDS - about to go to RightHandSide RHS (sweep E=%e RHO=%e)\n",sweep.rhs[2][ENG],sweep.rhs[2][RHO]);
		  RightHandSide (&sweep, Dts, nbeg, nend, dt, grid,1);
		  printf ("In UDS - back from RHS RightHandSide RHS sweep E=%e RHO=%e\n",sweep.rhs[2][ENG],sweep.rhs[2][RHO]);
		  
	  }
	  else {RightHandSide (&sweep, Dts, nbeg, nend, dt, grid,0);
	  }
		  

	  if (i==2 || j==2) printf ("In UDS2h  %e j=%3i i=%3i current E %e \n",g_time,j,i,UU[0][2][2][ENG]);

    /* -- Update:  U = U + dt*R -- */

      for ((*ip) = nbeg; (*ip) <= nend; (*ip)++) { 
		  if (j==2 && i==2) {if (*ip==2) printf ("RHS energy=%e %i %i\n",sweep.rhs[*ip][ENG],j,i);}
        NVAR_LOOP(nv) UU[k][j][i][nv] += sweep.rhs[*ip][nv];
      }
	  
	  if (i==2 || j==2) printf ("In UDS2i  %e j=%3i i=%3i current E %e \n",g_time,j,i,UU[0][2][2][ENG]);
	  


    /* -- Compute inverse hyperbolic time step - */

      #if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
      if (g_intStage == 1){
        inv_dl = GetInverse_dl(grid);
        for ((*ip) = nbeg; (*ip) <= nend; (*ip)++) { 
          C_dt[k][j][i] += 0.5*(Dts->cmax[(*ip)-1] + Dts->cmax[*ip])*inv_dl[*ip];
        }
      }
      #else
      inv_dl = GetInverse_dl(grid);
      for ((*ip) = nbeg-1; (*ip) <= nend; (*ip)++) { 
        Dts->invDt_hyp = MAX(Dts->invDt_hyp, Dts->cmax[*ip]*inv_dl[*ip]);
      }
      #endif
	  
    }
  }


  printf ("In UDS2  %e current E %e \n",g_time,UU[0][2][2][ENG]);


  

/* -------------------------------------------------------------------
   8. Reduce dt for dimensionally unsplit schemes.
   ------------------------------------------------------------------- */

#if (DIMENSIONAL_SPLITTING == NO) && (DIMENSIONS > 1)
  
  if (g_intStage == 1){
    DOM_LOOP(k,j,i) Dts->invDt_hyp = MAX(Dts->invDt_hyp, C_dt[k][j][i]);
    Dts->invDt_hyp /= (double)DIMENSIONS;
  }
#endif
  
  
  
}
