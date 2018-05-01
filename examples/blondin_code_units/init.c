/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
	double temp,cisosqrd,cent_mass,rho_alpha,disk_mdot,rho_0,r_0,mu,tx,r;
	
	cent_mass=g_inputParam[CENT_MASS];  //Central mass
	rho_alpha=g_inputParam[RHO_ALPHA];  //Drop off exponent of density
	disk_mdot=g_inputParam[DISK_MDOT];  //disk accretion rate - needed for the temperature
	rho_0=g_inputParam[RHO_0];          //density at r_0 - typically set to R_IC
	r_0=g_inputParam[R_0];              //r_0 - typically set to R_IC
    tx=g_inputParam[T_x];               //temperature of xray source - needed here to compute compton temperature
	
	
	
    r=x1*UNIT_LENGTH; //We will do things in cgs here - then convert back
	
	
//we work out the temperaure in an accretion disk at r=x1	

  temp=pow(((3.0*cent_mass*disk_mdot*CONST_G)/(8.0*CONST_PI*CONST_sigma*pow((r*sin(x2)),3.0))),0.25);
/*  mu=MeanMolecularWeight(v); This is how it should be done - but we are ionized and get the wrong answer */
  
  mu=0.6;
  
  
//we need the isothermal sound speed to work out the hydrostatic disk structure 
	
//  cisosqrd=g_gamma*temp/(KELVIN*mu);
  
  cisosqrd=g_gamma*CONST_kB*temp/CONST_mp;
  
	  
  v[RHO]=rho_0*pow((r/r_0),-1.0*rho_alpha); //midplane density at this radius
  
  
  v[RHO]=v[RHO]*exp(-1.0*CONST_G*cent_mass*pow((r*cos(x2)),2)/(cisosqrd*2.0*pow(r,3))); //hydrostatic strucuture
 
 
 
  if (v[RHO]<1e-22)  //Set a lower density bound
  {
  	v[RHO]=1e-22; 
  	temp=tx/4.0;    //If we have set the lower density - set the teperature to the compton temperature
}

printf ("r=%e v[rho]=%e %e\n",r,v[RHO],v[RHO]/UNIT_DENSITY);


	v[RHO]=v[RHO]/UNIT_DENSITY; //Scale to code units


  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3]=sqrt((CONST_G*cent_mass*sin(x2)*sin(x2))/r)/UNIT_VELOCITY;

  #if HAVE_ENERGY
//  v[PRS] = cisosqrd*v[RHO]/g_gamma;
  v[PRS]=v[RHO]*temp/(KELVIN*mu);   //This converts temperature (in K) to pressure in code units.
  #endif
  v[TRC] = 0.0;
  


  #if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;

  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double rho_0,r_0,rho_alpha,cent_mass;
  
//Get some parameters from the input file
    
rho_0=g_inputParam[RHO_0];
r_0=g_inputParam[R_0];
rho_alpha=g_inputParam[RHO_ALPHA];
cent_mass=g_inputParam[CENT_MASS];



  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (side == 0) 
  {    /* -- check solution inside domain -- */
	TOT_LOOP(k,j,i)
	{
		if (j==NX2+1)  //This should be the last 'real' theta bin - before the ghost zones.
		{
			d->Vc[RHO][k][j][i]=rho_0*pow((x1[i]/r_0),-1.0*rho_alpha);  //Set density at the midplane
			d->Vc[VX1][k][j][i]=0.0;									//Set radial velocity at the midplane to zero
		    d->Vc[VX3][k][j][i]=sqrt((CONST_G*cent_mass*sin(x2[j])*sin(x2[j]))/x1[i]);	  //Set v_phi to keplarian		
		}
		if (d->Vc[RHO][k][j][i] < 1.e-22) //Set a lower density thruoghout the domain
		{
	         d->Vc[RHO][k][j][i] = 1.e-22;
		 } 
	 }	
  }
  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -This is the midplane- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
  

}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
	double cent_mass;
	
	cent_mass=g_inputParam[CENT_MASS];
	
	
#if GEOMETRY == CARTESIAN
	return -1.0*cent_mass*CONST_G/sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY == CYLINDRICAL 
	return -1.0*cent_mass*CONST_G/sqrt(x1*x1 + x2*x2);
#elif GEOMETRY == SPHERICAL 
	return -1.0*cent_mass*CONST_G/x1;
#endif
}
#endif
