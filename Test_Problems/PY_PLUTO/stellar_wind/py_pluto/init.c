

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
	
    double temp,ciso2,cent_mass,rho_alpha,disk_mdot,rho_0,r_0,r,v_0;
    double H0; //The density scale weight in an isothermal wind
    double g0; //base of the wind
    
    cent_mass=g_inputParam[CENT_MASS];
	rho_0=g_inputParam[RHO_0];
	r_0=g_inputParam[R_0];
    v_0=g_inputParam[V_0];
    
    temp=g_inputParam[T_ISO];
    
    ciso2=(CONST_Rgas*temp/0.6);
    
    
    r=x1*UNIT_LENGTH;

        
    temp=1e5;    
    g0=CONST_G*cent_mass/r_0/r_0;    
    H0=CONST_Rgas*temp/0.6/g0;
        
    v[RHO]=rho_0*exp(-1*(r-r_0)/H0*r_0/r)/UNIT_DENSITY;

//    v[RHO]=rho_0*exp(-1.0*CONST_G*cent_mass/ciso2*((1./r_0)-(1./r)))/UNIT_DENSITY;
    
    v[VX1]=v_0*r/r_0/UNIT_VELOCITY;
    v[VX2]=0.0;
    v[VX3]=0.0;
    
    if (v[RHO]*UNIT_DENSITY<1e-30)
    {
        v[RHO]=1e-30/UNIT_DENSITY;
    }

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
    
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    double rho_0,r_0,r;
    double v_0;
    double v_l,v_r,dvdr,cent_mass;
    double ciso2,temp,massloss;
    
    rho_0=g_inputParam[RHO_0];
    v_0=g_inputParam[V_0];
    r_0=g_inputParam[R_0];
    cent_mass=g_inputParam[CENT_MASS];
    temp=g_inputParam[T_ISO];
    

    
    ciso2=(CONST_Rgas*temp/0.6);
    
    /* -- check solution inside domain -- */
    
 /*   
   	TOT_LOOP(k,j,i) //Now compute dvdr for use in line driving calculations
    {  
    
    
        printf ("BOOM1 %i %e %e %e\n",i,dvdr_array[i],d->Vc[RHO][0][0][i],d->Vc[VX1][0][0][i]);
    
}*/
    
    
        
  	DOM_LOOP(k,j,i)
  	{
        if (i==IBEG)  //First cell in the r-dimension
        {            
            d->Vc[RHO][k][j][i]=rho_0/UNIT_DENSITY; //Set the density BC
        }
  		if (d->Vc[RHO][k][j][i]*UNIT_DENSITY < 1.e-30) //Set a lower density throughout the domain
  		{
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
        if (i==0) //The first 'ghost' cell
        {
            v_l=(d->Vc[VX1][k][j][1]-d->Vc[VX1][k][j][0])/(grid->x[IDIR][1]-grid->x[IDIR][0]);
            v_l=d->Vc[VX1][k][j][0]-v_l*(grid->x[IDIR][0]-grid->xl[IDIR][0]);
            v_r=(d->Vc[VX1][k][j][1]-d->Vc[VX1][k][j][0])/(grid->x[IDIR][1]-grid->x[IDIR][0]);
            v_r=d->Vc[VX1][k][j][0]+v_r*(grid->xr[IDIR][0]-grid->x[IDIR][0]);
        }
        else if (i==IEND+2)
        {
            v_l=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            v_l=d->Vc[VX1][k][j][i]-v_l*(grid->x[IDIR][i]-grid->xl[IDIR][i]); 
            v_r=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            v_r=d->Vc[VX1][k][j][i]+v_r*(grid->xr[IDIR][i]-grid->x[IDIR][i]);                     
        }
        else
        {
            v_r=(d->Vc[VX1][k][j][i+1]-d->Vc[VX1][k][j][i])/(grid->x[IDIR][i+1]-grid->x[IDIR][i]);
            v_r=d->Vc[VX1][k][j][i]+v_r*(grid->xr[IDIR][i]-grid->x[IDIR][i]);            
            v_l=(d->Vc[VX1][k][j][i]-d->Vc[VX1][k][j][i-1])/(grid->x[IDIR][i]-grid->x[IDIR][i-1]);
            v_l=d->Vc[VX1][k][j][i]-v_l*(grid->x[IDIR][i]-grid->xl[IDIR][i]);          
        } 
        dvdr_array[i]=(v_r-v_l)/(grid->xr[IDIR][i]-grid->xl[IDIR][i]);	
    }    
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
  double k_rad,alpha_rad,***M;
  double test,M_test;
  double sigma_fc,mu_fc,D_fc,v_r,Rstar;
  
  
  
  k_rad=0.381;
  alpha_rad=0.595;
  
  M = GetUserVar("M");
  
    M_max=g_inputParam[M_rad];
    
    r=x1*UNIT_LENGTH;
    v_r=v[VX1]*UNIT_VELOCITY;
     
    rho=v[RHO]*UNIT_DENSITY;     
    L_x=g_inputParam[Lum_x]; 
    L_UV=g_inputParam[Lum_uv];
    
    Rstar=g_inputParam[R_0];
    
    F_UV=L_UV/4.0/CONST_PI/r/r;
    F_x=L_x/4.0/CONST_PI/r/r;
    
#if PY_CONNECT
    v_th=sqrt(2.*CONST_kB*py_temp[k][j][i]/CONST_amu);        
#else
    v_th=sqrt(2.*CONST_kB*g_inputParam[T_ISO]/CONST_amu);        
#endif
    
    sigma_e=CONST_sigmaT/CONST_amu/1.18;
    
    t=sigma_e*rho*v_th/fabs(dvdr_array[i]*UNIT_VELOCITY/UNIT_LENGTH); 
      
    M[k][j][i]=k_rad*pow(t,-1.0*alpha_rad);
    if (M[k][j][i]>M_max) M[k][j][i]=M_max;
    
    sigma_fc=fabs(r/v_r*(dvdr_array[i])*UNIT_VELOCITY/UNIT_LENGTH-1.);
    mu_fc=sqrt(1-(Rstar/r)*(Rstar/r));    
    D_fc=pow((1.+sigma_fc),(alpha_rad+1))-pow((1+sigma_fc*mu_fc*mu_fc),(alpha_rad+1));
    D_fc=D_fc/(1-mu_fc*mu_fc);
    D_fc=D_fc/(alpha_rad+1.)/sigma_fc;
    D_fc=D_fc/(pow((1+sigma_fc),alpha_rad));
     M[k][j][i]= M[k][j][i]*D_fc;
#if PY_CONNECT
    if (g_rad[0][k][j][i]==12345678 && g_rad[1][k][j][i]==12345678 && g_rad[2][k][j][i]==12345678) //First time thruogh, we will be using an approximation
    {
        g[IDIR]=(1.+M[k][j][i])*sigma_e*F_UV/CONST_c/UNIT_ACCELERATION;
        g[JDIR] = 0.0;
        g[KDIR] = 0.0;
    }
    else
    {
        g[IDIR]=g_rad[0][k][j][i];
        g[JDIR]=g_rad[1][k][j][i];
        g[KDIR]=g_rad[2][k][j][i];  
        M[k][j][i]=g[IDIR]/(sigma_e*F_UV/CONST_c/UNIT_ACCELERATION)-1.;     
    }
#else
    g[IDIR]=(1+M[k][j][i])*sigma_e*F_UV/CONST_c/UNIT_ACCELERATION;
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
	
    cent_mass=g_inputParam[CENT_MASS];
    
//       printf ("BOOM %e %e %e\n",cent_mass,CONST_G,(cent_mass*CONST_G/69570000000.0/69570000000.0));
    g_const=CONST_G/(UNIT_ACCELERATION*UNIT_LENGTH*UNIT_LENGTH/UNIT_MASS);
//    printf ("POTENTIAL %e\n",x1);
    
	cent_mass=g_inputParam[CENT_MASS]/UNIT_MASS;
//    printf ("BLAH g_isoSoundSpeed=%e\n",g_isoSoundSpeed) ; 
	
//    printf ("BOOM %e\n",(cent_mass*CONST_G));
    
#if GEOMETRY == CARTESIAN
	return -1.0*cent_mass*CONST_G/sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY == CYLINDRICAL 
	return -1.0*cent_mass*CONST_G/sqrt(x1*x1 + x2*x2);
#elif GEOMETRY == SPHERICAL
	return -1.0*cent_mass*g_const/x1;
    
//	return -1.0/x1;
    
#endif
}
#endif

