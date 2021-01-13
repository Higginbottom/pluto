

double average_dt;

#include "pluto.h"
/* ///////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////// */




/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
	
    double temp,ciso2,cent_mass,rho_alpha,disk_mdot,rho_0,r_0,r,v_0,z,z_0;
    double H0; //The density scale weight in an isothermal wind
    double g0; //base of the wind
    double theta; //angle to the cell - used to compute components of g field
    
    double g_r,g_z;
    
    cent_mass=g_inputParam[CENT_MASS];
	rho_0=g_inputParam[RHO_0];
    v_0=g_inputParam[V_0];
    z_0=g_inputParam[Z_0];
    
    temp=g_inputParam[T_ISO];
    
    ciso2=(CONST_Rgas*temp/0.6);
    g_isoSoundSpeed=sqrt(ciso2)/UNIT_VELOCITY;
    
    r=x1*UNIT_LENGTH;
    z=x2*UNIT_LENGTH;


    theta=atan(z/r);

    g_r=CONST_G*cent_mass/(r*r+z*z)*cos(theta);
    g_z=CONST_G*cent_mass/(r*r+z*z)*cos(theta);
        
//    temp=1.32e8;    
    g0=CONST_G*cent_mass/r_0/r_0;    
    H0=CONST_Rgas*temp/0.6/g0;
//    v[RHO]=rho_0*exp(-1*(r-r_0)/H0*r_0/r)/UNIT_DENSITY;

//    v[RHO]=rho_0*exp(-1.0*CONST_G*cent_mass/ciso2*((1./r_0)-(1./r)))/UNIT_DENSITY;
    
    v[VX1]=0.0;
    v[VX2]=v_0*z/z_0/UNIT_VELOCITY;
    v[VX3]=0.0;
   
//    printf ("BOOM %e %e\n",z,v[VX2]);
//    printf ("BOOM %e %e %e %e %e\n",r,z,g_r,g_z,v[VX3]);
    
//    if (v[RHO]*UNIT_DENSITY<1e-30)
//    {
//        v[RHO]=1e-30/UNIT_DENSITY;
//    }

//    v[VX3]=v[VX3]=sqrt(1.0/x1);
    
//    printf ("%e %e %e\n",r/UNIT_LENGTH,r,v[RHO]*UNIT_DENSITY);
//    printf ("BLAH g_isoSoundSpeed=%e\n",g_isoSoundSpeed) ; 


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
    double rho_0,v_0,cent_mass;
    int i,j,k;
    double theta;
    double g_r,g_z,r,z,temp,ciso2;
    
	rho_0=g_inputParam[RHO_0];
    v_0=g_inputParam[V_0];
    cent_mass=g_inputParam[CENT_MASS];
    
    temp=g_inputParam[T_ISO];
    ciso2=(CONST_Rgas*temp/0.6);
    
    printf ("CISO=%e\n",sqrt(ciso2));
    
    double *r_grid = grid->x[IDIR]; //The centre of the grid
    double *z_grid = grid->x[JDIR]; //The centre of the grid
    
    double *z_glob = grid->x_glob[JDIR];
    
    
  	DOM_LOOP(k,j,i)
  	{
        z=z_grid[j]*UNIT_LENGTH;
        r=r_grid[i]*UNIT_LENGTH;
        
        theta=atan(z/r);
        
        
        g_r=CONST_G*cent_mass/(r*r+z*z)*cos(theta);
        g_z=CONST_G*cent_mass/(r*r+z*z)*sin(theta);
        d->Vc[VX3][k][j][i]=sqrt(g_r*r)/UNIT_VELOCITY;
        
        
        
        
        if (j==JBEG && z_grid[j]==z_glob[2])  //First cell in the rdimension
            
            
//        	if (j==grid->np_int[JDIR]+1 && 	(fabs(x2_glob[grid->np_int_glob[JDIR]+1]-x2[j])/x2_glob[grid->np_int_glob[JDIR]+1]) < 1e-20)  //This should be the last 'real' theta bin - before the ghost zones.
            
            
        {              
            d->Vc[RHO][k][j][i]=rho_0/UNIT_DENSITY; //Set the density BC
            d->Vc[VX2][k][j][i]=v_0/UNIT_VELOCITY; //v_z            
        }
        else
        {        
            d->Vc[RHO][k][j][i]=rho_0*exp(-1*CONST_G*cent_mass*z*z/ciso2/2/r/r/r)/UNIT_DENSITY; //Isothermal atmpsphere
                
//        printf ("BOOM1 %i %i %i z=%e r=%e,g_r=%e,g_z=%e %e\n",k,j,i,z,r,g_r,g_z,exp(-1*CONST_G*cent_mass*z*z/ciso2/2/r/r/r));
                
                
//        if (d->Vc[RHO][k][j][i] < 1.e-20/UNIT_DENSITY)
            d->Vc[RHO][k][j][i] = 1.e-30/UNIT_DENSITY;				
      	}

        
  	 }
     
     
     
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *
 *********************************************************************** */
{
    int   i, j, k, nv;
    double *x1 = grid->x[IDIR]; //The centre of the grid
    double *xr1 = grid->xr[IDIR]; //right hand side (larger than centre)
    double *xl1 = grid->xl[IDIR]; //Left hand side (smaller than centre)

    
    double *r_grid = grid->x[IDIR]; //The centre of the grid
    double *z_grid = grid->x[JDIR]; //The centre of the grid
    
    double *z_glob = grid->x_glob[JDIR];
    

    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    double rho_0,r_0;
    double v_0;
    double vr_l,vr_r,dvdr,cent_mass;
    double vz_l,vz_u;
    
    double ciso2,temp,massloss;
    
    double z,r;
    
    double g_r,theta;
    
    rho_0=g_inputParam[RHO_0];
    v_0=g_inputParam[V_0];
//    r_0=g_inputParam[R_0];
    cent_mass=g_inputParam[CENT_MASS];
    temp=g_inputParam[T_ISO];
    

    
    ciso2=(CONST_Rgas*temp/0.6);
    
    
    
    //Set BCs in the computational domain - mainly in the first z
            
  	DOM_LOOP(k,j,i)
  	{
        
 //       print ("TEST %e %e %e %e\n",z_grid[j],z_glob[0],z_glob[1],z_glob[2]);
        
        if (j==JBEG && z_grid[j]==z_glob[2])  //First cell in the z-dimension (taking account of split domain...)
        {           
            print ("BOOM %e\n",z_grid[j]);
             
            z=z_grid[j]*UNIT_LENGTH;
            r=r_grid[i]*UNIT_LENGTH;
        
            theta=atan(z/r);
        
        
            g_r=CONST_G*cent_mass/(r*r+z*z)*cos(theta);
            
           
            d->Vc[VX3][k][j][i]=sqrt(g_r*r)/UNIT_VELOCITY; //Set phi velocity at the base to keplarian
            
            d->Vc[VX1][k][j][i]=0.0; //no radial velocity
            d->Vc[RHO][k][j][i]=rho_0/UNIT_DENSITY; //set density to BC
            
        }
  		if (d->Vc[RHO][k][j][i]*UNIT_DENSITY < 1.e-30) //Set a lower density throughout the domain
  		{
//            printf ("resetting %i %i %i %e\n",k,j,i, d->Vc[RHO][k][j][i]*UNIT_DENSITY);
  	         d->Vc[RHO][k][j][i] = 1.e-30/UNIT_DENSITY;				
  		 }
  	 }

     
     if (side == X1_BEG) //boundary condition at the start of the radial zone
     {
        dvdr=(d->Vc[VX1][0][0][IBEG+1]-d->Vc[VX1][0][0][IBEG])/(x1[IBEG+1]-x1[IBEG]); //Calculate the grasdient of velocity in the first cell
        d->Vc[VX1][0][0][1]=d->Vc[VX1][0][0][IBEG]-dvdr*(x1[IBEG]-x1[1]); //Use that gradient to compute velocity in first two cells
        d->Vc[VX1][0][0][0]=d->Vc[VX1][0][0][IBEG]-dvdr*(x1[IBEG]-x1[0]);

        //Compute a hydrostatic atmosphere for the density of the first two ghost zones

        d->Vc[RHO][0][0][0]=rho_0*exp(-1.0*CONST_G*cent_mass/ciso2*((1./x1[IBEG]/UNIT_LENGTH)-(1./x1[0]/UNIT_LENGTH)))/UNIT_DENSITY;
        d->Vc[RHO][0][0][1]=rho_0*exp(-1.0*CONST_G*cent_mass/ciso2*((1./x1[IBEG]/UNIT_LENGTH)-(1./x1[1]/UNIT_LENGTH)))/UNIT_DENSITY;
        massloss=rho_0*d->Vc[VX1][0][0][IBEG];
        
        
         BOX_LOOP(box,k,j,i) //Apply the normal outflow BC for the other two compnents of velocity
         {
             d->Vc[VX2][k][j][i]=d->Vc[VX2][k][j][IBEG];             
             d->Vc[VX3][k][j][i]=d->Vc[VX3][k][j][IBEG];                                     
         }
     }
     
     if (side == X1_END) //boundary condition at the end of the radial zone
     {
         
         //as with the inner ghost zones - extrapolate velocities based on the last gradient
         
         dvdr=(d->Vc[VX1][0][0][IEND]-d->Vc[VX1][0][0][IEND-1])/(x1[IEND]-x1[IEND-1]);
         d->Vc[VX1][0][0][IEND+2]=d->Vc[VX1][0][0][IEND]+dvdr*(x1[IEND+2]-x1[IEND]);
         d->Vc[VX1][0][0][IEND+1]=d->Vc[VX1][0][0][IEND]+dvdr*(x1[IEND+1]-x1[IEND]); 
         
         //But this time extrapolate density to maintain mass loss rate.
         
          d->Vc[RHO][0][0][IEND+1]=d->Vc[RHO][0][0][IEND]*d->Vc[VX1][0][0][IEND]/d->Vc[VX1][0][0][IEND+1];
          d->Vc[RHO][0][0][IEND+2]=d->Vc[RHO][0][0][IEND]*d->Vc[VX1][0][0][IEND]/d->Vc[VX1][0][0][IEND+2];

          //check against negative densities

          if (d->Vc[RHO][0][0][IEND+1]<0.0) d->Vc[RHO][0][0][IEND+1]=d->Vc[RHO][0][0][IEND];
          if (d->Vc[RHO][0][0][IEND+2]<0.0) d->Vc[RHO][0][0][IEND+2]=d->Vc[RHO][0][0][IEND+1];

         BOX_LOOP(box,k,j,i) //apply normal BCs 
         {
             d->Vc[VX2][k][j][i]=d->Vc[VX2][k][j][IEND];             
             d->Vc[VX3][k][j][i]=d->Vc[VX3][k][j][IEND];                                     
         }
     }
  

  
     
   	TOT_LOOP(k,j,i) //Now compute dvdr for use in line driving calculations
    {
//        if (i==50) printf ("i %i j %i r %e zl %e z %e zr %e v_z %e",i,j,grid->x[IDIR][i],grid->xl[JDIR][j],grid->x[JDIR][j],grid->xr[JDIR][j],d->Vc[VX2][k][j][i]);

        if (i==0) //The first 'ghost' cell in the radial direction
        {
            vr_l=(d->Vc[VX1][k][j][1]-d->Vc[VX1][k][j][0])/(grid->x[IDIR][1]-grid->x[IDIR][0]);
            vr_l=d->Vc[VX1][k][j][0]-vr_l*(grid->x[IDIR][0]-grid->xl[IDIR][0]);
            vr_r=(d->Vc[VX1][k][j][1]-d->Vc[VX1][k][j][0])/(grid->x[IDIR][1]-grid->x[IDIR][0]);
            vr_r=d->Vc[VX1][k][j][0]+vr_r*(grid->xr[IDIR][0]-grid->x[IDIR][0]);
        }
        
        else if (j==0) //The first 'ghost' cell in the z direction direction
        {
            vz_l=(d->Vc[VX2][k][1][i]-d->Vc[VX2][k][0][i])/(grid->x[JDIR][1]-grid->x[JDIR][0]);
            vz_l=d->Vc[VX2][k][0][i]-vz_l*(grid->x[JDIR][0]-grid->xl[JDIR][0]);
            vz_u=(d->Vc[VX2][k][1][i]-d->Vc[VX2][k][0][i])/(grid->x[JDIR][1]-grid->x[JDIR][0]);
            vz_u=d->Vc[VX2][k][0][i]+vz_u*(grid->xr[JDIR][0]-grid->x[JDIR][0]);
        }
        
        
        else if (i==IEND+2)
        {
            vr_l=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            vr_l=d->Vc[VX1][k][j][i]-vr_l*(grid->x[IDIR][i]-grid->xl[IDIR][i]); 
            vr_r=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            vr_r=d->Vc[VX1][k][j][i]+vr_r*(grid->xr[IDIR][i]-grid->x[IDIR][i]);                     
        }
        
        else if (j==JEND+2) //Thr last 'ghost' cell in the z direction..
        {
            vz_l=(d->Vc[VX2][k][j][i]-d->Vc[VX2][k][j-1][i])/(grid->x[JDIR][j]-grid->x[JDIR][j-1]);
            vz_l=d->Vc[VX2][k][j][i]-vz_l*(grid->x[JDIR][j]-grid->xl[JDIR][j]); 
            vz_u=(d->Vc[VX2][k][j][i]-d->Vc[VX2][k][j-1][i])/(grid->x[JDIR][j]-grid->x[JDIR][j-1]);
            vz_u=d->Vc[VX2][k][j][i]+vz_u*(grid->xr[JDIR][j]-grid->x[JDIR][j]);                     
        }
        
        
        else
        {
            vr_r=(d->Vc[VX1][k][j][i+1]-d->Vc[VX1][k][j][i])/(grid->x[IDIR][i+1]-grid->x[IDIR][i]);
            vr_r=d->Vc[VX1][k][j][i]+vr_r*(grid->xr[IDIR][i]-grid->x[IDIR][i]);            
            vr_l=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            vr_l=d->Vc[VX1][k][j][i]-vr_l*(grid->x[IDIR][i]-grid->xl[IDIR][i]);  
            
            
            
            
            vz_u=(d->Vc[VX2][k][j+1][i]-d->Vc[VX2][k][j][i])/(grid->x[JDIR][j+1]-grid->x[JDIR][j]); //gradient of velocity across upper cell boundary
//            if (i==50) printf (" grad_upper %e ",vz_u);
 
            vz_u=d->Vc[VX2][k][j][i]+vz_u*(grid->xr[JDIR][j]-grid->x[JDIR][j]);   //vertical velocity at top of cell
                      
            vz_l=(d->Vc[VX2][k][j][i]-d->Vc[VX2][k][j-1][i])/(grid->x[JDIR][j]-grid->x[JDIR][j-1]); //gradient of velocity across lower cell boundary
//            if (i==50) printf (" grad_lower %e ",vz_l);
            
            vz_l=d->Vc[VX2][k][j][i]-vz_l*(grid->x[JDIR][j]-grid->xl[JDIR][j]); //vertical velocity at bottom of cell
            
                    
        } 
        dvr_dr_array[k][j][i]=(vr_r-vr_l)/(grid->xr[IDIR][i]-grid->xl[IDIR][i]);
        dvz_dz_array[k][j][i]=(vz_u-vz_l)/(grid->xr[JDIR][j]-grid->xl[JDIR][j]);	
//        if (i==50) printf (" vz_l %e vz_u %e dvz_dz %e\n",vz_l,vz_u,dvz_dz_array[k][j][i]);
        
//        printf ("BOOM %i %i %i %e %e %e %e %e %e\n",k,j,i,vz_u,vz_l,grid->xr[JDIR][j],grid->xl[JDIR][j],(grid->xr[JDIR][j]-grid->xl[JDIR][j]),dvz_dz_array[k][j][i]);
        	
    }
//   exit(0);    
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
  double A,cent_mass,Gamma,g_const,L_UV,L_x;
  double F_UV,F_x,temp,r;
  double A_test,g_r,sigma_e;
  double dvdr,t,rho,v_th,M_max;
  double k_r,alpha_r,***M;
  double test,M_test;
  double sigma_fc,mu_fc,D_fc,v_r,Rstar;
  double disk_mdot,flux,g_z,theta,z;
  
  
//  printf ("Arrived\n");
  r=x1*UNIT_LENGTH;
  z=x2*UNIT_LENGTH;
  
  theta=atan(z/r);
  
  M_max=g_inputParam[M_rad];
  
  M = GetUserVar("M");
  
  k_r=g_inputParam[k_rad];
  alpha_r=g_inputParam[alpha_rad];
  rho=v[RHO]*UNIT_DENSITY;     
  
  
  disk_mdot=g_inputParam[M_acc];
  cent_mass=g_inputParam[CENT_MASS];
  
#if PY_CONNECT
    v_th=sqrt(2.*CONST_kB*py_temp[k][j][i]/CONST_amu);        
#else
    v_th=sqrt(2.*CONST_kB*g_inputParam[T_ISO]/CONST_amu);        
#endif
  
  
  
  g_z=sin(theta)*CONST_G*cent_mass/(r*r+z*z);
  
  
  //Work out the temperature of the disk below this element
  temp=pow(((3.0*cent_mass*disk_mdot*CONST_G)/(8.0*CONST_PI*CONST_sigma*pow((r),3.0))),0.25);
  //and the flux from that disk element
  flux=CONST_sigma*pow(temp,4.);
  //this is the electron scattering opacity per unit mass
  sigma_e=CONST_sigmaT/CONST_amu/1.18;
  
  t=sigma_e*rho*v_th/fabs(dvz_dz_array[k][j][i]*UNIT_VELOCITY/UNIT_LENGTH); 
//  printf ("Here %e %e %e\n",t,v_th,1./fabs(dvz_dz_array[k][j][i]*UNIT_VELOCITY/UNIT_LENGTH));
    
  M[k][j][i]=k_rad*pow(t,-1.0*alpha_rad);
  if (M[k][j][i]>M_max) M[k][j][i]=M_max;
  
//  printf ("BOOM %e\n",M[k][j][i]);
  
//  printf ("BOOM %e\n",(M_max*sigma_e*flux/CONST_c)/UNIT_ACCELERATION);
 
    g[IDIR]=0.0; //acceleration in the omega direction
    g[JDIR]=(M[k][j][i]*sigma_e*flux/CONST_c)/UNIT_ACCELERATION; //and in the z direction
    g[KDIR]=0.0; //and in the phi direction

//   printf ("%e %e %e\n",z,g[JDIR]*UNIT_ACCELERATION,g_z);


//        printf ("BOOM3 %e %e %e %e %e\n",r,temp,flux,g[JDIR]*UNIT_ACCELERATION,g_z);
    
    
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
	
    cent_mass=g_inputParam[CENT_MASS];
    
//       printf ("BOOM %e %e %e\n",cent_mass,CONST_G,(cent_mass*CONST_G/69570000000.0/69570000000.0));
    g_const=CONST_G/(UNIT_ACCELERATION*UNIT_LENGTH*UNIT_LENGTH/UNIT_MASS);
//    printf ("POTENTIAL %e\n",x1);
    
	cent_mass=g_inputParam[CENT_MASS]/UNIT_MASS;
//    printf ("BLAH g_isoSoundSpeed=%e\n",g_isoSoundSpeed) ; 
	
//    printf ("BOOM %e\n",(cent_mass*CONST_G));
    
#if GEOMETRY == CARTESIAN
	return -1.0*cent_mass*g_const/sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY == CYLINDRICAL     
	return -1.0*cent_mass*g_const/sqrt(x1*x1 + x2*x2);
#elif GEOMETRY == SPHERICAL
	return -1.0*cent_mass*g_const/x1;
//	return -1.0/x1;
    
#endif
}
#endif

