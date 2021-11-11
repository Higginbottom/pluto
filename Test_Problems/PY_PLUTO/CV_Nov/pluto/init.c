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
	double disk_trunc_rad,T_iso;
    double x,z;
    double v_x,v_z;
	
	cent_mass=g_inputParam[CENT_MASS];  //Central mass
	rho_alpha=g_inputParam[RHO_ALPHA];  //Drop off exponent of density
	disk_mdot=g_inputParam[DISK_MDOT];  //disk accretion rate - needed for the temperature
	rho_0=g_inputParam[RHO_0];          //density at r_0 - typically set to R_IC
	r_0=g_inputParam[R_0];              //r_0 - typically set to R_IC
    tx=g_inputParam[T_x];               //temperature of xray source - needed here to compute compton temperature
    T_iso=g_inputParam[T_ISO];               //temperature of xray source - needed here to compute compton temperature
//    printf ("BOOM T_iso %e\n",T_iso);
	
    r=x1*UNIT_LENGTH; //We will do things in cgs here - then convert back
	
    
    
//we work out the temperaure in an accretion disk at r=x1	

#if EOS != ISOTHERMAL
  temp=pow(((3.0*cent_mass*disk_mdot*CONST_G)/(8.0*CONST_PI*CONST_sigma*pow((r*sin(x2)),3.0))),0.25);
#else
  temp=T_iso;
  
  g_isoSoundSpeed=sqrt(CONST_Rgas*temp/0.6)/UNIT_VELOCITY;
  
  g_isoSoundSpeed=14.*1000.*100./UNIT_VELOCITY; //14 km/s the value from PSD98
  
  
#endif
/*  mu=MeanMolecularWeight(v); This is how it should be done - but we are ionized and get the wrong answer */
 
  mu=0.6;
  
  dvds_setup_flag=0;
  
//we need the isothermal sound speed to work out the hydrostatic disk structure 
	
//  cisosqrd=g_gamma*temp/(KELVIN*mu);
  
  cisosqrd=g_isoSoundSpeed*g_isoSoundSpeed*UNIT_VELOCITY*UNIT_VELOCITY;
//  printf ("BOOM %e %e\n",temp,sqrt(cisosqrd));
 
//      v[RHO]=rho_0*pow((r/r_0),-1.0*rho_alpha); //midplane density at this radius
//      v[RHO]=v[RHO]*exp(-1.0*CONST_G*cent_mass*pow((r*cos(x2)),2)/(cisosqrd*2.0*pow(r,3))); //hydrostatic strucuture

  v[RHO]=rho_0*exp(-1*CONST_G*cent_mass/2./cisosqrd/r/tan(x2)/tan(x2)); //The expression from PSD98

 
  if (v[RHO]<1e-20)  //Set a lower density bound
  {
  	v[RHO]=1e-20; 
  	temp=tx/4.0;    //If we have set the lower density - set the teperature to the compton temperature
}


	v[RHO]=v[RHO]/UNIT_DENSITY; //Scale to code units

    x=x1*sin(x2)*UNIT_LENGTH;
    z=x1*cos(x2)*UNIT_LENGTH;

 //   v_x=x/1e2/UNIT_VELOCITY;
//    v_z=0.0;
        
        //z/1e2/UNIT_VELOCITY;
        
       
	v[VX1]=v_x*sin(x2)+v_z*cos(x2);
	v[VX2]=v_x*cos(x2)-v_z*sin(x2); 
        


//  v[VX1] = r/1e2/UNIT_VELOCITY;
  v[VX1] = 0.0;
   v[VX2] = 0.0;
 //  v[VX1] = r/1e2/UNIT_VELOCITY;
   
//  v[VX2] = x2*1e5/UNIT_VELOCITY;
   v[VX3]=sqrt((CONST_G*cent_mass*sin(x2)*sin(x2))/r)/UNIT_VELOCITY;
// v[VX3]=sqrt(CONST_G*cent_mass/r)/sin(x2)/UNIT_VELOCITY;
 
 
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
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
    double r,theta;
    int i,j,k;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    
    DOM_LOOP(k,j,i){
		r = x1[i]*UNIT_LENGTH;
		theta= x2[j];
//        printf ("r= %e theta=%e Density=%e iso_soundspeed=%e\n",r,theta,d->Vc[RHO][k][j][i]*UNIT_DENSITY,g_isoSoundSpeed*UNIT_VELOCITY);


}      
    
    
    
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
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
	  
  double *x2_glob = grid->x_glob[JDIR];
	  
  double  r, theta;
  double rho_0,r_0,rho_alpha,cent_mass,disk_trunc_rad,tx;
  double v_new;
//Get some parameters from the input file
    
rho_0=g_inputParam[RHO_0];
r_0=g_inputParam[R_0];
rho_alpha=g_inputParam[RHO_ALPHA];
cent_mass=g_inputParam[CENT_MASS];
tx=g_inputParam[T_x];               //temperature of xray source - needed here to compute compton temperature


  
  if (side == 0) 
  {    /* -- check solution inside domain -- */
	DOM_LOOP(k,j,i)
	{
	if (j==grid->np_int[JDIR]+1 && 	(fabs(x2_glob[grid->np_int_glob[JDIR]+1]-x2[j])/x2_glob[grid->np_int_glob[JDIR]+1]) < 1e-20)  //This should be the last 'real' theta bin - before the ghost zones.
		{
			r = x1[i]*UNIT_LENGTH;
			theta= x2[j];
            
            
            d->Vc[VX2][k][j][i]=d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]/((rho_0*pow((r/r_0),-1.0*rho_alpha))/UNIT_DENSITY); //conserve momentum
//            printf ("%e %e\n",d->Vc[RHO][k][j][i],(rho_0*pow((r/r_0),-1.0*rho_alpha))/UNIT_DENSITY);
            d->Vc[RHO][k][j][i]=(rho_0*pow((r/r_0),-1.0*rho_alpha))/UNIT_DENSITY;  //Set density at the midplane
			d->Vc[VX1][k][j][i]=0.0;									//Set radial velocity at the midplane to zero
		    d->Vc[VX3][k][j][i]=(sqrt((CONST_G*cent_mass*sin(theta)*sin(theta))/r))/UNIT_VELOCITY;	  //Set v_phi to keplarian in code units			 
		}
//    	if (i==2) 
//    		{
//                if (d->Vc[RHO][k][j][i]*UNIT_DENSITY<1e-14)
//                {
//                    d->Vc[RHO][k][j][i]=1e-14/UNIT_DENSITY;  //Set density at the midplane
//                }
	 
//    		}    
		if (d->Vc[RHO][k][j][i]*UNIT_DENSITY < 1.e-20) //Set a lower density throughout the domain
		{   
//            printf ("BOOM %e",d->Vc[VX1][k][j][i]);
            d->Vc[VX1][k][j][i]=d->Vc[VX1][k][j][i]/(1.e-20/UNIT_DENSITY/d->Vc[RHO][k][j][i]); //conserve momentum
            d->Vc[VX2][k][j][i]=d->Vc[VX2][k][j][i]/(1.e-20/UNIT_DENSITY/d->Vc[RHO][k][j][i]); //conserve momentum
//:wq
            d->Vc[VX3][k][j][i]=d->Vc[VX3][k][j][i]/(1.e-20/UNIT_DENSITY/d->Vc[RHO][k][j][i]); //conserve momentum           
//            printf (" %e\n",d->Vc[VX1][k][j][i]);
            
	         d->Vc[RHO][k][j][i] = 1.e-20/UNIT_DENSITY;	
             
             
             
//             d->Vc[PRS][k][j][i]=d->Vc[RHO][k][j][i]*tx/(KELVIN*0.6);			
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


void VGradCalc(const Data *d, Grid *grid) 
{
    int i,j,k,itest,jtest;
    double vx1_l,vx1_r,vx2_l,vx2_u;
    double v11[2],v12[2],v22[2],v21[2];
    double v1,v2,dx1,dx2;
    double v11_f[2],v12_f[2],v22_f[2],v21_f[2];
    double unit_flux[2];
    
    double x11[2],x22[2],flux_x,flux_z;
    double loc[2],ans1[2],ans2[2],mod_flux;
    
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    
    double *x1l = grid->xl[IDIR];
    double *x2l = grid->xl[JDIR];
    double *x3l = grid->xl[KDIR];
    
    double *x1r = grid->xr[IDIR];
    double *x2r = grid->xr[JDIR];
    double *x3r = grid->xr[KDIR];
    
    double x,z;
    double vx1,vz1,vx2,vz2;
    
    double ds,maxds,temp;
    
    int flag,iangle;
        
    if (dvds_setup_flag==0) //This is the first time here, so we need to set up a few things that wont change
	{
		printf ("Setting up dvds calculation arrays\n");
		//Firstly, we are able to set up the various arrays to hold the offsets for each angular bin		
    	dvds_r_offset = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
    	dvds_t_offset = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
    	dvds_mod_offset = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);	
		
		//And we can also set up the array for the dvds array - one value for each flux bin
    	dvds_array = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);	
					
		for (iangle=0;iangle<NFLUX_ANGLES;iangle++)
		{
	       	DOM_LOOP(k,j,i) //Now compute dvdr for use in line driving calculations
	    	{        
	        	//The points at which interpolations are carried out never change
        
	        	x11_interp[0][k][j][i]=(x1[i-1]+x1[i])/2.0*UNIT_LENGTH;
	        	x11_interp[1][k][j][i]=(x2[j-1]+x2[j])/2.0;
        
	        	x22_interp[0][k][j][i]=(x1[i+1]+x1[i])/2.0*UNIT_LENGTH;
	        	x22_interp[1][k][j][i]=(x2[j+1]+x2[j])/2.0;
        
	        	//We can also compute the second location to compute the velocity in each cell - that just depends on the flux
	        	//First we work out how far we want to move in this cell
                
	        	maxds=fabs(x22_interp[0][k][j][i]-x11_interp[0][k][j][i]); //The x1 distance between intepolation centres in real units.
       
	        	#if GEOMETRY == SPHERICAL //Approximate linear distance across the interpolation centres
	        		if (maxds>(x1[i]*UNIT_LENGTH*fabs(x22_interp[1][k][j][i]-x11_interp[1][k][j][i])))
	        		{
	            		maxds=x1[i]*UNIT_LENGTH*fabs(x22_interp[1][k][j][i]-x11_interp[1][k][j][i]);
	        		}
	        	#endif        
        
	        		maxds/=10.; //This should mean than the distance we move to compute dv/ds is within the interpolation grid
        
	        	//We now need the direction to move for the disk flux - that depends on the flux - we first work out the modulus of  the flux
				mod_flux=sqrt(pow(flux_x_UV[iangle][k][j][i],2)+pow(flux_y_UV[iangle][k][j][i],2)+pow(flux_z_UV[iangle][k][j][i],2));
			
	        	if (mod_flux==0.0) //There is no flux in this cell - capture the issue with a flag
	        	{
	            	dvds_mod_offset[iangle][k][j][i]=-999; //the fact that this is negative acts as the flag
	        	}               
	        	else
	        	{
	            	flux_x=flux_x_UV[iangle][k][j][i]/mod_flux;
	            	flux_z=flux_z_UV[iangle][k][j][i]/mod_flux;        
	            	#if GEOMETRY == SPHERICAL
	                x=x1[i]*sin(x2[j])*UNIT_LENGTH; //Get the cell centre in cartesian coords
	                z=x1[i]*cos(x2[j])*UNIT_LENGTH;    
	                dx1=maxds*flux_x;  //This is a move in the x direction
	                dx2=maxds*flux_z; //This is a move in the z direction
	                ds=sqrt(dx1*dx1+dx2*dx2); //This is belt and braces - it *should* be equal to maxds                 
	         	   //Go back to spherical coordinates because the interpolation routine works in r,theta and store the locations on a global array for constant use      
	                dvds_r_offset[iangle][k][j][i]=sqrt((x+dx1)*(x+dx1)+(z+dx2)*(z+dx2));
	                dvds_t_offset[iangle][k][j][i]=atan((x+dx1)/(z+dx2));
	                dvds_mod_offset[iangle][k][j][i]=ds;
	            	#endif	
	        	}
	    	} //End of DOM loop for initialising the first time round - everything after this point is done every time
		}
		printf ("Finished creating the dvds arrays\n");
	    dvds_setup_flag=1;
    }
    
//    printf ("Making velocity gradients\n");

	for (iangle=0;iangle<NFLUX_ANGLES;iangle++)
    {
	   	DOM_LOOP(k,j,i) //Now compute dvdr for use in line driving calculations
	    {         
	        //Firstly, compute the velocity at the interpolation vertices vertices
	        //If we simply compute the mean velocity at four points equidistant from the surrounding cell centers, thats fine - the locations dont matter.
        
	        v11[0]=(d->Vc[VX1][k][j-1][i-1]+d->Vc[VX1][k][j-1][i]+d->Vc[VX1][k][j][i-1]+d->Vc[VX1][k][j][i])/4.0;
	        v11[1]=(d->Vc[VX2][k][j-1][i-1]+d->Vc[VX2][k][j-1][i]+d->Vc[VX2][k][j][i-1]+d->Vc[VX2][k][j][i])/4.0;

	        v12[0]=(d->Vc[VX1][k][j+1][i-1]+d->Vc[VX1][k][j][i-1]+d->Vc[VX1][k][j+1][i]+d->Vc[VX1][k][j][i])/4.0;
	        v12[1]=(d->Vc[VX2][k][j+1][i-1]+d->Vc[VX2][k][j][i-1]+d->Vc[VX2][k][j+1][i]+d->Vc[VX2][k][j][i])/4.0;
        
	        v22[0]=(d->Vc[VX1][k][j+1][i]+d->Vc[VX1][k][j+1][i+1]+d->Vc[VX1][k][j][i+1]+d->Vc[VX1][k][j][i])/4.0;
	        v22[1]=(d->Vc[VX2][k][j+1][i]+d->Vc[VX2][k][j+1][i+1]+d->Vc[VX2][k][j][i+1]+d->Vc[VX2][k][j][i])/4.0;
        
	        v21[0]=(d->Vc[VX1][k][j][i+1]+d->Vc[VX1][k][j-1][i+1]+d->Vc[VX1][k][j-1][i]+d->Vc[VX1][k][j][i])/4.0;
	        v21[1]=(d->Vc[VX2][k][j][i+1]+d->Vc[VX2][k][j-1][i+1]+d->Vc[VX2][k][j-1][i]+d->Vc[VX2][k][j][i])/4.0;
        
        
	        //Now get the locations of the interpolation vertices - these are stored at the start
        
	        x11[0]=x11_interp[0][k][j][i];
	        x11[1]=x11_interp[1][k][j][i];
	        x22[0]=x22_interp[0][k][j][i];
	        x22[1]=x22_interp[1][k][j][i];

			// We will be computing the velocity at two points in the direction of the flux - we need to work out a sensible step across the interpolation space
	        //We first get the velocity at the cell centre - this *should* match the actual cell center velocity, so we have a test
    
	        loc[0]=x1[i]*UNIT_LENGTH;
	        loc[1]=x2[j];
    
	        bilinear(x11,x22,v11,v12,v21,v22,loc,ans1,0);
                    
            
	        //Get the velocities at first location in cartesian coordinates    
	        vx1=(ans1[0]*UNIT_VELOCITY*sin(x2[j])+ans1[1]*UNIT_VELOCITY*cos(x2[j]));
	        vz1=(ans1[0]*UNIT_VELOCITY*cos(x2[j])-ans1[1]*UNIT_VELOCITY*sin(x2[j]));            
        
			//Extract the second location which is stored for each angular bin
			
	        loc[0]=dvds_r_offset[iangle][k][j][i];
	        loc[1]=dvds_t_offset[iangle][k][j][i];
	        ds=dvds_mod_offset[iangle][k][j][i];
        
	        if (ds==-999) //This is a flag that there is no flux in this cell - so dvds cant be calculated
	        {
	            dvds_array[iangle][k][j][i]=-999; //A flag to later routines that we dont have a flux here, so the gradient is undefined
	        }
	        else
	        {           
                
	    		//Compute the velocity at the new offset point                
        
	            bilinear(x11,x22,v11,v12,v21,v22,loc,ans2,0);
        
	        	//Get the velocities at second location in x,z coordinates    
	            vx2=(ans2[0]*UNIT_VELOCITY*sin(loc[1])+ans2[1]*UNIT_VELOCITY*cos(loc[1]));
	            vz2=(ans2[0]*UNIT_VELOCITY*cos(loc[1])-ans2[1]*UNIT_VELOCITY*sin(loc[1])); 
                
	            //we now have the two velocities, seperated by a distance ds along the direction of of flux. We now need the dot products of each along the flux direction.
        
				mod_flux=sqrt(pow(flux_x_UV[iangle][k][j][i],2)+pow(flux_y_UV[iangle][k][j][i],2)+pow(flux_z_UV[iangle][k][j][i],2));
		
            	flux_x=flux_x_UV[iangle][k][j][i]/mod_flux;
            	flux_z=flux_z_UV[iangle][k][j][i]/mod_flux; 
		
   	            v1=flux_x*vx1+flux_z*vz1;  
	            v2=flux_x*vx2+flux_z*vz2;  

	            dvds_array[iangle][k][j][i]=fabs((v2-v1)/ds); //These are the velocity gradients in physical quantities.

	        }
//			if (iangle==0 && k==0 && i==2 && j==2)
//			{
//	            printf ("TEST               iang=%i k=%i j=%i i=%i dvds_array %e flux %e %e %e %e ds %e\n",iangle,k,j,i,dvds_array[iangle][k][j][i],flux_x_UV[iangle][k][j][i],flux_y_UV[iangle][k][j][i],flux_z_UV[iangle][k][j][i],mod_flux,ds);
				
//			}
	        if (dvds_array[iangle][k][j][i]!=dvds_array[iangle][k][j][i])
	        {
	            printf ("BOOM we have a nan in cycle %li iang=%i k=%i j=%i i=%i dvds_array %e flux %e %e %e %e ds %e\n",g_stepNumber,iangle,k,j,i,dvds_array[iangle][k][j][i],flux_x_UV[iangle][k][j][i],flux_y_UV[iangle][k][j][i],flux_z_UV[iangle][k][j][i],mod_flux,ds);       
	            exit(0);
	        } 
	    } //This the end of the DOM loop
	} //This is the end of the loop over angular bins
}
    


void bilinear (double x11[2],double x22[2],double v11[2],double v12[2],double v21[2],double v22[2],double test[2],double ans[2],int flag)
{
    double fracx1,fracx2;
    double temp1,temp2;
    

    
    
    fracx1 = (test[0] - x11[0]) / (x22[0] -x11[0]);  //x1 direction
    fracx2 = (test[1] - x11[1]) / (x22[1] -x11[1]);  //x2 direction
    if (flag==1)
    
    printf ("BOOM frac1=%e frac2=%e\n",fracx1,fracx2);
    
    temp1= (1. - fracx1) * v11[0] + fracx1 * v21[0]; //vx1 interpolation at the bottom of the cell
    temp2= (1. - fracx1) * v12[0] + fracx1 * v22[0]; //vx1 interpolation at the bottom of the cell
    
    if (flag==1)
    
    printf ("BOOM temp1=%e temp2=%e\n",temp1,temp2);
    
    ans[0]=(1. - fracx2) * temp1 + fracx2 * temp2;  //vx1 interpolation in the theta direction
    
    temp1= (1. - fracx1) * v11[1] + fracx1 * v21[1]; //vx2 interpolation at the bottom of the cell
    temp2= (1. - fracx1) * v12[1] + fracx1 * v22[1]; //vx2 interpolation at the top of the cell
    
    ans[1]=(1. - fracx2) * temp1 + fracx2 * temp2;  //vx2 interpolation in the theta direction
    if (flag==1)
    printf ("v11 %e v12 %e v21 %e v22 %e\n",v11[0],v12[0],v21[0],v22[0]);
    
}





#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3,int i,int j,int k)
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
	double rad, rs,theta;
	double Lx,Luv,flux;
	double rho,nH,ne,g_r,g_t,g_p;
	double g_x,g_y,g_z;
    double sigma_e;
    double t_UV;
    double v_th;
    double krad_UV,alpha_UV;
    double M_UV;
	double M_UV_array[MPOINTS];
    double eta_max,tau_max,M_max,T_iso;
	int iangle,ii;
    
    T_iso=g_inputParam[T_ISO];               //isothermal temperature
    krad_UV=k_UV_array[k][j][i];               
    alpha_UV=alpha_UV_array[k][j][i];                
//    printf ("BOOM %i\n",MPOINTS );


	for (ii=0;ii<MPOINTS;ii++)
	{
		M_UV_array[ii]=M_UV_fit[ii][0][j][i];
	}

	
	

	sigma_e=CONST_sigmaT/CONST_amu/1.18;
//    v_th=4.2e5;
    rho = v[RHO]*UNIT_DENSITY;  //density of the current cell in physical units
//    krad=0.2;
//   alpha=0.6;
    M_max=4400.;  
    
    
    v_th=  pow ((2. * CONST_kB * T_iso / CONST_mp), 0.5);     //We need the thermal velocity for hydrogen
//    printf ("BOOM %e %e\n",T_iso,v_th);
//   eta_max=pow((M_max/(krad_UV*(1.-alpha_UV))),(1./alpha_UV));
    
//    printf ("eta_max=%e\n",eta_max);
    	
	Lx=g_inputParam[L_star]*g_inputParam[f_x];  //central lum in cgs
	Luv=g_inputParam[L_star]*g_inputParam[f_uv];  //central lum in cgs
	  	
	nH=rho/(1.43*CONST_mp);   //Work out hydrogen number density assuming stellar abundances
	ne=1.21*nH;             //electron number density assuming full ionization
	
	/* first we get the radial distance from the central source*/
	
	#if GEOMETRY == CARTESIAN
		rs = sqrt(x1*x1 + x2*x2 + x3*x3); /* spherical radius in cart. coords */
	#elif GEOMETRY == CYLINDRICAL
		rs = sqrt(x1*x1 + x2*x2); /* spherical radius in cyl. coords */
	#elif GEOMETRY == SPHERICAL
		rs = x1*UNIT_LENGTH; /* spherical radius in sph. coords */
	#endif
		
	/* and also the angle of the cell */ 
		
	#if GEOMETRY == CARTESIAN
		theta = atan(x1/x2); /* spherical radius in cart. coords */
	#elif GEOMETRY == CYLINDRICAL
		theta = atan(x1/x2); /* spherical radius in cyl. coords */
	#elif GEOMETRY == SPHERICAL
		theta = x2; /* spherical radius in sph. coords */
	#endif	
		
	//copy across the cartesian accelerations from python
	
    flux=Luv/4.0/CONST_PI/rs/rs;
    
#if PY_CONNECT
    #if PY_RAD_DRIV==ACCELERATIONS
        if (g_rad[0][k][j][i]==12345678 && g_rad[1][k][j][i]==12345678 && g_rad[2][k][j][i]==12345678) //First time thruogh, we will be using an approximation
        {
            g[IDIR]=(2000.)*sigma_e*flux/CONST_c/UNIT_ACCELERATION;
            g[JDIR] = 0.0;
            g[KDIR] = 0.0;
        }
        else
        {
        	g_x=g_rad[0][k][j][i];
        	g_y=g_rad[1][k][j][i];
        	g_z=g_rad[2][k][j][i];
    		g[IDIR] = g_x*sin(theta)+g_z*cos(theta); 
    		g[JDIR] = g_x*cos(theta)-g_z*sin(theta);
    		g[KDIR] = g_y;	
    //        printf ("BOOM %e %e\n",g_x*sin(theta)+g_z*cos(theta),(2000.)*sigma_e*flux/CONST_c/UNIT_ACCELERATION);  
        }
    #endif
    #if PY_RAD_DRIV==FLUXES
		//The total acceleration is the sum over all directional bins
		g[IDIR]=0.0;
		g[JDIR]=0.0;
		g[KDIR]=0.0;
		for (iangle=0;iangle<NFLUX_ANGLES;iangle++)
		{				
//			printf ("Doing angle %i flux %e  dvds %e\n",iangle,flux_r_UV[iangle][k][j][i],dvds_array[iangle][k][j][i]);
//			printf ("Doing angle %i flux %e  dvds %e\n",iangle,flux_t_UV[iangle][k][j][i],dvds_array[iangle][k][j][i]);
			
			//We need to compute the force multiplier, using the velocity gradient in this cell.
		    if (dvds_array[iangle][k][j][i]>0.0)
		    {
//				printf ("interp1\n");
		        t_UV=   sigma_e * rho * v_th / dvds_array[iangle][k][j][i];  
//				printf ("interp2\n");
				
				M_UV=linterp(log10(t_UV),t_fit,M_UV_array,MPOINTS);
				
				
//				printf ("BOOM %f j=%i i=%i iangle=%i\n",M_UV,j,i,iangle);
				
//		        M_UV=krad_UV*pow(t_UV,(alpha_UV));
		        if (M_UV>M_max)
		            M_UV=M_max;
		    }
		    else
		    {				
		        M_UV=M_max; //This is a bit of a fudge - the only way this should be negative is if there is no flux..
		    }		
	        if (g_stepNumber>1)
	        {
	            g[IDIR]=g[IDIR]+(1+M_UV)*sigma_e*flux_r_UV[iangle][k][j][i]/CONST_c/UNIT_ACCELERATION;
	            g[JDIR]=g[JDIR]+(1+M_UV)*sigma_e*flux_t_UV[iangle][k][j][i]/CONST_c/UNIT_ACCELERATION;
	        }
	        else
	        {
	            g[IDIR]=g[IDIR]+(1.)*sigma_e*flux_r_UV[iangle][k][j][i]/CONST_c/UNIT_ACCELERATION;
	            g[JDIR]=g[JDIR]+(1.)*sigma_e*flux_t_UV[iangle][k][j][i]/CONST_c/UNIT_ACCELERATION;
	        }
	        g[KDIR] = 0.0; 
//			printf ("Finished doing angle %i flux %e  dvds %e\n",iangle,flux_r_UV[iangle][k][j][i],dvds_array[iangle][k][j][i]);
		}
    #endif    
    
#else
    g[IDIR]=(2000.)*sigma_e*flux/CONST_c/UNIT_ACCELERATION;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;    
#endif    
    
        
		

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
	double cent_mass,g_const;
	
	
	cent_mass=g_inputParam[CENT_MASS]/UNIT_MASS;  //central mass in code units
	g_const=CONST_G/(UNIT_ACCELERATION*UNIT_LENGTH*UNIT_LENGTH/UNIT_MASS);  //grav constant in code units
		
	
#if GEOMETRY == CARTESIAN
	return -1.0*cent_mass*g_const/sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY == CYLINDRICAL 
	return -1.0*cent_mass*g_const/sqrt(x1*x1 + x2*x2);
#elif GEOMETRY == SPHERICAL 
	return -1.0*cent_mass*g_const/x1;
#endif
}
#endif

double linterp (double x, double xarray[], double yarray[], int nelem)
{

	int ii;
	double ans,dydx;
	
	ii=0;
	while (xarray[ii]<x)
	{
		ii++;
	}
//	printf ("%e is between %e and %e\n",x,xarray[ii-1],xarray[ii]);
	
	if (yarray[ii-1]<-1)
	{
//		printf ("We are off the end of the curve M=%e\n",pow(10,yarray[ii-1]));
		ans=0.0;
	}
	else if (ii>nelem)
	{
		ans=pow(10,yarray[ii-1]);
	}		
	else
	{
		dydx=(yarray[ii]-yarray[ii-1])/(xarray[ii]-xarray[ii-1]);
//	printf ("ylo=%e yhi=%e dydx=%e \n",pow(10,yarray[ii-1]),pow(10,yarray[ii]),dydx);
		ans=pow(10,yarray[ii-1]+(x-xarray[ii-1])*dydx);
	}
	
	if (ans!=ans)
	printf ("BOOM NAN %e xlo %e xhi %e ylo %e yhi %e\n",ans,xarray[ii-1],xarray[ii],yarray[ii-1],yarray[ii]);

	
	return(ans);
}




