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

#define LINELENGTH 800



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
	double gr,gt,gp;
	#endif
	

	if (flag==0) //Initialising first time round
    {
    	printf ("Entering py_heatcool - first time round\n");        
    	DOM_LOOP(k,j,i){
        
#if COOLING == BLONDIN
			d->comp_h_pre[k][j][i]=1.0;
			d->comp_c_pre[k][j][i]=1.0;
			d->xray_h_pre[k][j][i]=1.0;
			d->line_c_pre[k][j][i]=1.0;
			d->brem_c_pre[k][j][i]=1.0;
#endif
#if PY_RAD_DRIV == ACCELERATIONS         
			g_rad[0][k][j][i]=1e-99;
			g_rad[1][k][j][i]=1e-99;
			g_rad[2][k][j][i]=1e-99;
#endif	
		}
		printf ("Initialised\n");
	}	
	else //We will be reading in some updated numbers
	{
		printf ("Reading in updated heating and cooling prefactors\n");
		fptr = fopen ("prefactors.dat", "r");
		printf ("Opened prefactor file\n");
		if (fptr==NULL)
		{
		    DOM_LOOP(k,j,i){
#if COOLING == BLONDIN                
				d->comp_h_pre[k][j][i]=1.0;
				d->comp_c_pre[k][j][i]=1.0;
				d->xray_h_pre[k][j][i]=1.0;
				d->line_c_pre[k][j][i]=1.0;
				d->brem_c_pre[k][j][i]=1.0;
#endif                
    			g_rad[0][k][j][i]=1e-99;
    			g_rad[1][k][j][i]=1e-99;
    			g_rad[2][k][j][i]=1e-99;		
			}
			printf ("NO prefactor file\n");
		}
		else
		{
	    DOM_LOOP(k,j,i){ //Initialise
#if COOLING == BLONDIN                            
			d->comp_h_pre[k][j][i]=-1.0;
			d->comp_c_pre[k][j][i]=-1.0;
			d->xray_h_pre[k][j][i]=-1.0;
			d->line_c_pre[k][j][i]=-1.0;
			d->brem_c_pre[k][j][i]=-1.0;
#endif   
#if PY_RAD_DRIV == ACCELERATIONS         
			g_rad[0][k][j][i]=1e-99;
			g_rad[1][k][j][i]=1e-99;
			g_rad[2][k][j][i]=1e-99;
#endif
		    }
		fgets (aline, LINELENGTH, fptr);
		icount=0;
		while (fgets (aline, LINELENGTH, fptr) != NULL)	
		{		
			if 
				#if (BODY_FORCE & VECTOR)
				((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le %le %le %le", &ii, &rcen, &jj, &thetacen, &dens,
				&comp_h_pre, &comp_c_pre, &xray_h_pre, &brem_c_pre, &line_c_pre, &gr, &gt, &gp)) == 13)
				#else
				((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le",             &ii, &rcen, &jj, &thetacen, &dens,
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
#if COOLING == BLONDIN                                                    
						d->comp_h_pre[k][j][i]=comp_h_pre;
						d->comp_c_pre[k][j][i]=comp_c_pre;
						d->xray_h_pre[k][j][i]=xray_h_pre;
						d->line_c_pre[k][j][i]=line_c_pre;
						d->brem_c_pre[k][j][i]=brem_c_pre;
#endif
						#if (BODY_FORCE & VECTOR & PY_RAD_DRIV == ACCELERATIONS)
							g_rad[0][k][j][i]=gr; //The acceleration as computed by python
							g_rad[1][k][j][i]=gt;
							g_rad[2][k][j][i]=gp;
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
#if COOLING == BLONDIN                                                            
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
#endif
        
    
	
		printf ("Read in %i prefectors\n",icount);
	}
}
}



/*This routine is used when we are running isothermal runs with driving*/
#if EOS==ISOTHERMAL && PY_CONNECT
/* ********************************************************************* */
void read_py_iso_temp (Data *d, Grid *grid,int flag)
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
	double dens,temp,nh,ne,t_opt,t_UV,t_Xray,r;
	int icount;
	#if (BODY_FORCE & VECTOR)
	double gx,gy,gz;
	#endif
    
    
    double r1=2e12/UNIT_LENGTH;
    double t0=50000;
    double t1=35000; 
    double dtdr; 
    
    dtdr=(t1-t0)/(r1-grid->x[IDIR][0]); 
    
    
    temp=g_inputParam[T_ISO];
	
    if (flag==1)
    {
        printf ("We are restarting so expecting to read in temperatures\n");
        fptr = fopen ("py_pcon_data.dat", "r");
        if (fptr==NULL)
        {
            printf ("NO pcon file\n");
    	}
    	fgets (aline, LINELENGTH, fptr);
    	icount=0;
        k=0;
    	while (fgets (aline, LINELENGTH, fptr) != NULL)	
    	{		
    	    if ((nwords = sscanf (aline, "%d %d %le %le %le %le %le %le %le", &ii, &jj, &temp, &dens, &nh, &ne, &t_opt, &t_UV, &t_Xray)) == 9)
    	    {
                if (fabs((dens-d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY)/dens) > 1e-6)
                    {
                        printf ("NOT SAME  %e %i\n",fabs((dens-d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY)/dens),k);
                    }
                    else
                    {
//                        printf ("Reading into %li %li\n",jj+JBEG,ii+IBEG);
                        py_temp[k][jj+JBEG][ii+IBEG]=temp;
//                        printf ("TEMP READ IN %li %e %e\n",ii+IBEG,grid->x[IDIR][ii+IBEG]*UNIT_LENGTH,py_temp[k][jj+JBEG][ii+IBEG]);
                        
                    }
            }
            else
            {
                printf ("Incorrectly formatted\n");
            }
    	}
        TOT_LOOP(k,j,i) //Fill in the ends of the temperature array
        {
            if (i<IBEG)
            {
                py_temp[k][j][i]=py_temp[k][j][IBEG];
            }
            else if (i>IEND)
            {
                py_temp[k][j][i]=py_temp[k][j][IEND];
            }
        }
    }
    else
    {
    	printf ("First time round so setting temperatures to fixed T\n");
        TOT_LOOP(k,j,i)
        {
            if (grid->x[IDIR][i] < r1)
            {
                py_temp[k][j][i]=t0+(grid->x[IDIR][i]-grid->x[IDIR][0])*dtdr;
            }
            else
            {
                py_temp[k][j][i]=t1;//For the initial run - we have no driving                	
            }
            py_temp[k][j][i]=temp;
//            printf ("%i %e %e\n",i,grid->x[IDIR][i]*UNIT_LENGTH,py_temp[k][j][i]);            
        } 
               	
    }
    
    printf ("%li %li\n",IBEG,IEND);
    
//    TOT_LOOP(k,j,i)
//    {
//        printf ("We will be using ISO temperatures %i %e %e\n",i,grid->x[IDIR][i]*UNIT_LENGTH,py_temp[k][j][i]);       
//    }
    
    	
}

#endif

/*This routine is used when we are running isothermal runs with driving*/

/* ********************************************************************* */
void read_py_rad_driv (Data *d, Grid *grid,int flag)
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
	double rcen,thetacen,dens,gx,gy,gz;
	int icount;

	printf ("Reading in accelerations from python\n");
    if (flag==1)
    {
        fptr = fopen ("py_accelerations.dat", "r");
        if (fptr==NULL)
        {
            printf ("NO accelerations file\n");
    	}
    	fgets (aline, LINELENGTH, fptr);
    	icount=0;
    	while (fgets (aline, LINELENGTH, fptr) != NULL)	
    	{		
    	    if ((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le",  &ii, &rcen, &jj, &thetacen, &dens, &gx, &gy, &gz)) == 8)
    	    {
                if (fabs((dens-d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY)/dens) > 1e-6)
                    printf ("NOT SAME  %e %e %e %e\n",thetacen,dens,gx,d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY);
                else
                {
                    icount++;
    				g_rad[0][k][jj+JBEG][ii+IBEG]=gx/UNIT_ACCELERATION;
    				g_rad[1][k][jj+JBEG][ii+IBEG]=gy/UNIT_ACCELERATION;
    				g_rad[2][k][jj+JBEG][ii+IBEG]=gz/UNIT_ACCELERATION;
                }
                
            }
            else
            {
                printf ("Incorrectly formatted\n");
            }
    	}
        printf ("matched %i accelerations\n",icount);
    }
    else
    {
        DOM_LOOP(k,j,i){
        	printf ("WE ARE HERE - first time round\n");
    		g_rad[0][k][j][i]=12345678.; //For the initial run - we have no driving - set a flag
    		g_rad[1][k][j][i]=12345678.;
    		g_rad[2][k][j][i]=12345678.;		
    	}        
    }	
}



void read_py_fluxes (Data *d, Grid *grid)
/*
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
	int i,j;
	long int k;
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	long int ii,jj,nwords;
	double x1in,x2in;
    double kradin,alpharadin,rhoin;
    double krad,alpharad;
    double opt[4],UV[4],Xray[4];
	double temp;
	int iflux;
	int icount,iaxis;
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    double tol=1e-6;
    k=0;
    krad=g_inputParam[KRAD];               
    alpharad=g_inputParam[ALPHARAD]; 
	
	
	
	printf ("Arrived in read_py_fluxes IBEG=%li JBEG=%li\n",IBEG,JBEG);
	printf ("Reading in fluxes from python\n");

	for (iaxis=0;iaxis<3;iaxis++) //We loop over three axes - x,y,z in seperate files
	{
		if (iaxis==0) fptr = fopen ("directional_flux_x.dat", "r");
		else if (iaxis==1) fptr = fopen ("directional_flux_y.dat", "r");
		else if (iaxis==2) fptr = fopen ("directional_flux_z.dat", "r");
	    if (fptr==NULL)
	    {
	        printf ("No flux file\n");
		}
		fgets (aline, LINELENGTH, fptr);
		nwords = sscanf (aline, "%*s %ld ",  &ii);
		if (nwords == 1)
		{
			if (iaxis==0)
			{			
				NFLUX_ANGLES=ii;
				printf ("We have %i angular bins for flux\n",NFLUX_ANGLES);	
				printf ("Reading in X fluxes\n");
		    	flux_x_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data	
			}
			else if (iaxis==1)
			{
				if (ii!=NFLUX_ANGLES)
				{
					printf ("Y-flux file doesnt agree in NFLUX_ANGLES CRASH!!\n");
					exit(0);
				}
				else
				{	
			    	flux_y_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data
					printf ("Reading in Y fluxes\n");

				}			
			}
			else if (iaxis==2)
			{
				if (ii!=NFLUX_ANGLES)
				{
					printf ("Z-flux file doesnt agree in NFLUX_ANGLES CRASH!!\n");
					exit(0);
				}
				else
				{
			    	flux_z_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data
					printf ("Reading in Z fluxes\n");
					
				}			
			}
		}
		else
		{
			printf ("Flux header improperly formatted\n");
		}
		icount=0;
		while (fscanf(fptr,"%ld ",&ii) !=EOF) //Read the first element of a line, and check for EOF
		{
			fscanf(fptr,"%ld %le %le",&jj,&x1in,  &x2in); //Get the rest of the geometry stuff for this line
            if (fabs(1.-(x1in/UNIT_LENGTH/x1[ii+IBEG]))<tol && fabs(1.-(x2in/x2[jj+IBEG]))<tol)
			{
				for (iflux=0;iflux<NFLUX_ANGLES;iflux++)
				{
					fscanf (fptr,"%le",&temp); //read the flux in one angular bin at a time
					if (iaxis==0)
						flux_x_UV[iflux][k][jj+JBEG][ii+IBEG]=temp;
					if (iaxis==1)
						flux_y_UV[iflux][k][jj+JBEG][ii+IBEG]=temp;
					if (iaxis==2)
						flux_z_UV[iflux][k][jj+JBEG][ii+IBEG]=temp;
				}
			icount++;	
			}
			else
	        {
	        	printf ("missed tolerance\n");
	        }
		}
		printf ("Read %i fluxes for %i cells\n",NFLUX_ANGLES,icount);
	    fclose(fptr);
	}	
	
	
	flux_r_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	flux_t_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	flux_p_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	
//Make fluxes in r,theta,phi directions as well	
	
	for (iflux=0;iflux<NFLUX_ANGLES;iflux++)
	{ 
		DOM_LOOP(k,j,i)
		{
			flux_r_UV[iflux][k][j][i]=flux_x_UV[iflux][k][j][i]*sin(x2[j])+flux_z_UV[iflux][k][j][i]*cos(x2[j]);
			flux_t_UV[iflux][k][j][i]=flux_x_UV[iflux][k][j][i]*cos(x2[j])-flux_z_UV[iflux][k][j][i]*sin(x2[j]);
			flux_p_UV[iflux][k][j][i]=flux_y_UV[iflux][k][j][i];		
		}
	}
	
	
	k=0; //Need to reset this number to zero after the dom loop!
	
	
    if (krad==999 && alpharad==999)
    {
		printf ("Importing a force multiplier fit file\n");
		fptr=fopen ("M_UV_data.dat", "r");
	    if (fptr==NULL)
	    {
	        printf ("No force multiplier file\n");
		}
		fgets (aline, LINELENGTH, fptr);
		nwords = sscanf (aline, "%*s %ld ",  &ii);
		MPOINTS=ii;
		printf ("There are %i points in the force multiplier fits\n",MPOINTS);
		M_UV_fit=ARRAY_4D(MPOINTS,NX3_TOT, NX2_TOT, NX1_TOT, double);
		
		fscanf(fptr,"%*s "); //Get the rest of the geometry stuff for this line
		t_fit=calloc(MPOINTS, sizeof(double));
		
        for (iflux=0;iflux<MPOINTS;iflux++)
		{
			fscanf (fptr,"%le",&temp); //read the values of t for which the force multiplier is tabulates
			t_fit[iflux]=log10(temp);		
//			printf ("BOOM %e\n",t_fit[iflux]);	

		}
		
		while (fscanf(fptr,"%ld ",&ii) !=EOF) //Read the first element of a line, and check for EOF
		{
			fscanf(fptr,"%ld ",&jj); //Get the rest of the geometry stuff for this line
//			printf ("Processing cell %i %i\n",ii,jj);
            for (iflux=0;iflux<MPOINTS;iflux++)
			{
				fscanf (fptr,"%le",&temp); //read the flux in one angular bin at a time
//				printf ("INPUT %i %e\n",iflux,temp);
				M_UV_fit[iflux][k][jj+JBEG][ii+IBEG]=log10(temp);
//				M_UV_fit[iflux][k][jj+JBEG][ii+IBEG]=log10(0.59*pow(pow(10,t_fit[iflux]),-0.6));
//				if (jj+JBEG==2 && ii+IBEG==2)
//				{
//				printf ("BOOM %li %li %li %i %i %e\n",k,jj+JBEG,ii+IBEG,iflux,MPOINTS,M_UV_fit[iflux][k][jj+JBEG][ii+IBEG]);	
//				}
			}
			icount++;	
		}	    
		fclose(fptr);
    }
    else
    {
        printf ("Using k=%e alpha=%e for all cells\n",krad,alpharad);
    }
	/*        printf("We are importing k alpha\n");
        fptr = fopen ("k_alpha_fit.dat", "r");
        
    	fgets (aline, LINELENGTH, fptr);
    	icount=0;
    	while (fgets (aline, LINELENGTH, fptr) != NULL)	
    	{		
    	    if ((nwords = sscanf (aline, "%d %d %le %le %le %le ",  &ii, &jj, &x1in,  &x2in,  &kradin, &alpharadin)) == 6)
    	    {
    #if GEOMETRY == SPHERICAL
                if (fabs(1.-(x1in/UNIT_LENGTH/x1[ii+IBEG]))<tol && fabs(1.-(x2in/x2[jj+IBEG]))<tol)
                {
    //                    printf ("BOOM %e %e %d %d (%ld %ld)%e (%e) %e (%e)\n",1.-(rcen/UNIT_LENGTH/x1[ii+IBEG]), 1.-(thetacen/x2[jj+IBEG]),ii,jj,jj+JBEG,ii+IBEG,rcen,rcen/x1[ii+IBEG],thetacen,thetacen/x2[jj+JBEG]);
                    icount++;                                
                    k_UV_array[k][jj+JBEG][ii+IBEG]=kradin;
                    alpha_UV_array[k][jj+JBEG][ii+IBEG]=alpharadin;                                    
                }
                else
                {
                    printf ("missed tolerance\n");
                }
    #else
                printf ("k-alpha import in this geometry not currently supported")
    #endif                    
            }
            else
            {
                printf ("Incorrectly formatted %d words\n",nwords);
            }
    	}
        printf ("matched %i pairs of k-alpha\n",icount);         */

    
}




