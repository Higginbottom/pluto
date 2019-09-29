#include "cak.h"

//#define f_uv_lower 7.5e14
//#define f_uv_upper 3e16

#define f_uv_lower 7.5e14
#define f_uv_upper 3e16


double scale=1.0;

/* Main - top level routine. */
int main(int argc, char** argv)
{

  int nelem, nion, nlvl, ntran;
  FILE *line_list, *atom_models, *output, *spew;
  int i,j;
  int dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8;
  float dum10, dum11, dum12, dum13, dum14, dum15;
  double partition;
  double line_nu;
  double A_ul, B_ul, B_lu;
  double kappa_l, delta_doppler, flux_factor, test;
  int upper, lower;
  int icell,iband,itrans,iion,ilev;
  int n_used_lines;
  int v_uv_x; //An index used to decide if we are incremeneting the vis(0), uv(1)or xray(2) band 
  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel
  int my_nmin,my_nmax;
  int np_mpi_global,num_mpi_cells,num_mpi_extra,rank_global,ndo;
  double *M_array_transmit;
  double *M_array_transmit2;


#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif
np_mpi_global = np_mpi;
rank_global=my_rank;
  printf ("BLAH doing it myrank=%i/%i\n",my_rank,np_mpi);
  g_usedMemory=0.0;
  


  read_ionfs();  //Read in the ion fractions and also the 
    
  read_cont();  //Read in the continuum
  
  read_line_data(); //Read in all the line data

//  exp_lookup_init();

//	M_array=ARRAY_2D(ncells,nbands,double); This is the way we need to do it once we have bands!
	M_array=ARRAY_2D(ncells,3,double); //Three bands, vis, UV and Xray
	
	part=malloc(nions*sizeof(double));
	  
	printf ("Set up working arrays - currently using %f Mb\n",g_usedMemory/1e6);
	printf ("Calculating\n");
	
	printf ("icell=%i itrans=%i\n",ncells,nline_tot);
    
    my_nmin = 0;
    my_nmax = ncells;
#ifdef MPI_ON
  num_mpi_cells = floor (ncells / np_mpi_global);
  num_mpi_extra = ncells - (np_mpi_global * num_mpi_cells);
  
  printf ("BLAH num_mpi_cells=%i num_mpi_extra=%i\n",num_mpi_cells,num_mpi_extra);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
  if (rank_global < num_mpi_extra)
  {
    my_nmin = rank_global * (num_mpi_cells + 1);
    my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
  }
  else
  {
    my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra) * (num_mpi_cells);
    my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra + 1) * (num_mpi_cells);
  }
  ndo = my_nmax - my_nmin;
#endif
    
    
    
    
    
    
  printf ("This thread doing cells %i to %i\n",my_nmin,my_nmax);
    
//	for (icell=0; icell<ncells; icell++)  //We are gping to loop over all cells
    for (icell = my_nmin; icell < my_nmax; icell++)
	{
		n_used_lines=0;
/* first thing we need to do is set up the level populations for this cell*/
		for (iion=0;iion<nions;iion++)  //loop over all ions
		{
			part[iion]=0.0;
			for (ilev=0;ilev<ion_info_nlevels[iion];ilev++) //loop over all the levels for this ion
			{
				lev_pops[iion][ilev]=lev_weight[iion][ilev]*exp1(-1. * lev_energy[iion][ilev] *EV/KB/T_e[icell]); //Level population - based on boltzmann distribution
				part[iion] += lev_pops[iion][ilev]; 					
			}
			part[iion] = 1./part[iion];	
			for (ilev=0;ilev<ion_info_nlevels[iion];ilev++) //loop over all the levels for this ion
			{
				lev_pops[iion][ilev]*=part[iion]*ion_fracs[icell][iion]; //Level population - based on boltzmann distribution
			}			
		}					
/*We now have level populations and partition function for all the ions in this cell*/ 		
//		for (iband=0; iband<nbands; iband++) //qLoop over all the bands	
			for (itrans=0;itrans<nline_tot;itrans++) //Loop over all lines
			{
				if (trans_freq[itrans] > 0.0) 
				{
					iion = trans_ion[itrans]; //This gets the link into the various ion related arrays
					if (ion_fracs[icell][iion]>0.0)
					{
						if (trans_freq[itrans]<f_uv_lower) v_uv_x=0;
						else if (trans_freq[itrans]>f_uv_upper) v_uv_x=2;
						else v_uv_x=1;
						n_used_lines++;
						A_ul = trans_A_ul[itrans];
						B_ul = trans_B_ul[itrans];
						B_lu = trans_B_lu[itrans];
						lower= trans_lower[itrans];
						upper= trans_upper[itrans];
						line_nu= trans_freq[itrans];
						kappa_l = (B_lu*lev_pops[iion][lower] - B_ul*lev_pops[iion][upper]) * param1[i]; 
						delta_doppler = line_nu * v_thermal[icell]/C;
						if (J[icell]==0.0)
							flux_factor=0.0;
						else
							flux_factor = model_jnu(line_nu, trans_lfreq[itrans],icell,trans_iband[itrans])/J[icell];							
						if ( (test = kappa_l*t[icell][v_uv_x]/sigma_e[icell]) > 1.e-6) //If we have a large opacity
						{	
							M_array[icell][v_uv_x] += delta_doppler * flux_factor * (1. - exp1(-1.*test)) / t[icell][v_uv_x];  //increment the fore multiplier for this t
							if (!isfinite(M_array[icell][v_uv_x]))
							{
								printf("Non-finite M1. Abort.\n");
								printf("%i %i %i %g %g %g %g\n", icell,cell_i[icell],cell_j[icell],delta_doppler,flux_factor,exp1(-1.*test), t[icell][v_uv_x]);
								exit(0);
							}				
						}
						else
						{
							M_array[icell][v_uv_x] += delta_doppler * flux_factor * kappa_l / sigma_e[icell];						
							if (!isfinite(M_array[icell][v_uv_x]))
							{
								printf("Non-finite M2. Abort.\n");
								printf("%g %g %g\n", delta_doppler,flux_factor,kappa_l);
								exit(0);
							}
						}
					} //End of loop to only do computations if ion frac is greater than zero							
				} //End of loop to compute M if the trans freq is greater than zero 
				else
				{					
					printf ("Dodgy line %i ion %i lev %i-%i  freq %e ignoring\n",itrans,trans_ion[itrans],trans_upper[itrans],trans_lower[itrans],trans_freq[itrans]);
				}					
			}	//End of the loop over all lines for this ion 
			printf ("Cell %i %i lines used\n",icell,n_used_lines);
	}//End of the loop over all cells

	printf ("Done all the calculations\n");
    
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile
    MPI_Barrier (MPI_COMM_WORLD);
    
	M_array_transmit=malloc(ncells*3.*sizeof(double));
	M_array_transmit2=malloc(ncells*3.*sizeof(double));
    
    for (i=0;i<ncells;i++)
    {
        M_array_transmit[i]=M_array[i][0];
        M_array_transmit[i+ncells]=M_array[i][1];
        M_array_transmit[i+2*ncells]=M_array[i][2];
    }
    MPI_Barrier (MPI_COMM_WORLD);    
    MPI_Reduce (M_array_transmit, M_array_transmit2, ncells*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (M_array_transmit2, ncells*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
    for (i=0;i<ncells;i++)
    {
        M_array[i][0]=M_array_transmit2[i];
        M_array[i][1]=M_array_transmit2[i+ncells];
        M_array[i][2]=M_array_transmit2[i+2*ncells];
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
#endif
    
    
    
#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif


    if ((output = fopen("M_data.dat", "w")) == NULL)
      {
	printf("Cannot open M_data.dat\n");
	abort();
      }
	fprintf(output, "i j M_vis M_uv M_xray\n");
	  	  
	  
	for (icell = 0; icell < ncells; icell++)
	{
		fprintf(output, "%i %i %e %e %e", cell_i[icell], cell_j[icell], M_array[icell][0], M_array[icell][1], M_array[icell][2]);
		fprintf(output,"\n");
	}
    fclose(output);
#ifdef MPI_ON
  }
#endif
    
    
#ifdef MPI_ON
  MPI_Finalize ();
#endif

}


/**************************************************************************************/
//																					
// This bit of code reads in the ion details.	
//
// It populates the array ion_fracs[NELEMS][NIONS] which are the abundances of ions
//
/**************************************************************************************/


int read_ionfs()
{
  int nelem,i,j;
  FILE *abund,*pcon;
  char filename[100];
  char junk[LINELENGTH];
  float dum1,dum2,dum3,dum4,dum5,dum6,dum7;
  int idum1,idum2;
  char *found;
  
  
/* A loop over the elements - NB at present we need *all* the elements in place - in files which comtain the number of the element as the filename*/  


  if ((abund=fopen("py_ion_data.dat","r")) ==NULL)
  {
	  printf ("Cannot open py_ion_data.dat\n");
	  abort();
  }
  
	fscanf(abund, "%*s %g", &dum1);
	nions=dum1;
	printf ("We have %i ions\n",nions);
	
  	fgets(junk, 1000, abund); 
	
	
	ion_info_z=malloc(nions*sizeof(int));
	ion_info_state=malloc(nions*sizeof(int));
	ion_info_nlevels=malloc(nions*sizeof(int));
	ion_info_nlines=malloc(nions*sizeof(int));
	ion_info_xi=malloc(nions*sizeof(double));
	
	
	g_usedMemory+=4.*nions*sizeof(int);
	g_usedMemory+=nions*sizeof(double);

	
	for (i=0;i<nions;i++)
	{
		if (fgets (junk, LINELENGTH, abund) == NULL)
        {
         	printf ("Error reading ion data %i of %i\n",i,nions);
         	abort();
        }
		sscanf (junk, "%*s %*d %*s %d %d", &idum1, &idum2);
		ion_info_z[i]=idum1;
		ion_info_state[i]=idum2;
	}
	
	printf ("We have read in the ion data - last element is z=%i state=%i\n",ion_info_z[nions-1],ion_info_state[nions-1]);
		
	fscanf(abund, "%*s %d", &idum1);
	ncells=idum1;
	
	printf("We have %i cells\n",ncells);
	
	ion_fracs=ARRAY_2D(ncells,nions,double); //Set up the ion fracs array
		
	for (i=0;i<ncells;i++)
	{
		fscanf(abund, " %d", &idum1);
		fscanf(abund, " %d", &idum2);
		for (j=0;j<nions;j++)
		{
			fscanf(abund, "%e",&dum1);
			ion_fracs[i][j]=dum1;
		}
	}
    fclose(abund);  //Close the file
	
	printf ("Read ion data - currently using %f Mb\n",g_usedMemory/1e6);
	
	
	
    if ((pcon=fopen("py_pcon_data.dat","r")) ==NULL)
    {
  	  printf ("Cannot open py_pcon_data.dat\n");
  	  abort();
    }	
	
	if (fgets (junk, LINELENGTH, pcon) == NULL)
    {
     	printf ("Error reading first line of py_pcon_data.dat\n");
     	abort();
    }
		
	sscanf(junk, "%*s %d", &idum1);
	
	
	if (idum1!=ncells)
	{
		printf ("py_pcon_data has different number of cells - aborting");
		abort();
	}


	cell_i=malloc(ncells*sizeof(int)); 
	cell_j=malloc(ncells*sizeof(int)); 


	T_e=malloc(ncells*sizeof(double)); 
	nnh=malloc(ncells*sizeof(double));
	nne=malloc(ncells*sizeof(double));
	rho=malloc(ncells*sizeof(double));
	v_thermal=malloc(ncells*sizeof(double));
	sigma_e=malloc(ncells*sizeof(double));
	param1=malloc(ncells*sizeof(double));
	
	g_usedMemory+=8.*ncells*sizeof(double);
	g_usedMemory+=2.*ncells*sizeof(int);
	
	
	t=ARRAY_2D(ncells,3,double); //Three bands, vis, UV and Xray
	
	
	for (i=0;i<ncells;i++)
	{
		if (fgets (junk, LINELENGTH, pcon) == NULL)
        {
         	printf ("Error reading pcon data %i of %i\n",i,ncells);
         	abort();
        }
		sscanf (junk, "%d %d %e %e %e %e %e %e %e", &idum1,&idum2,&dum1, &dum2,&dum3,&dum4,&dum5,&dum6,&dum7);
		cell_i[i]=idum1;
		cell_j[i]=idum2;
		T_e[i]=dum1;
		rho[i]=dum2;
		nnh[i]=dum3;
		nne[i]=dum4;
		t[i][0]=dum5;	
		t[i][1]=dum6;	
		t[i][2]=dum7;	
		
		v_thermal[i] = pow((2. * KB * T_e[i]/MH),0.5);
		sigma_e[i] = SIGMA_T * nne[i] / rho[i]; //electron scattering sigma	
		param1[i]=HCLIGHTOVERFOURPI / rho[i] / v_thermal[i];
		
	}
	
	
	printf ("Read physical data - currently using %f Mb\n",g_usedMemory/1e6);
	
	
	return(0);
}
/***************************************************************************/
// Read_cont reads in the continuum flux in the cell. Originally, it was
// set up to read in from a cloudy file - 
/***************************************************************************/


int read_cont()
{  	
  FILE *contf, *output;
  char filename[100];
  char junk[1000];
  double dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9, dum10;
  int idum1,idum2;
  int i,j;
  
  
	if ((contf = fopen("py_spec_data.dat", "r")) == NULL)
	{
		printf("Cannot open py_spec_data.dat\n");
		abort();
	}
	
	fscanf(contf, "%*s %d", &idum1);
	nbands=idum1;
	
	printf ("We have %i bands\n",nbands);
	
	
	band_limits=malloc(nbands*sizeof(double));
	
	fscanf(contf, "%*s %d", &idum1);

	if (ncells!=idum1)
		{
			printf ("Big problem - number of cells in spec (%i) doesnt equal those in pcon (%i)\n",idum1,ncells);
			abort();
		}
	
	printf ("We have %i cells\n",ncells);
	
	for (j=0;j<nbands+1;j++)
	{
		fscanf(contf, "%le ",&dum1);
		band_limits[j]=dum1;
	}

	
	f1=ARRAY_2D(ncells,nbands,double);
	f2=ARRAY_2D(ncells,nbands,double);
	model=ARRAY_2D(ncells,nbands,int);	
	pl_w=ARRAY_2D(ncells,nbands,double);
	pl_alpha=ARRAY_2D(ncells,nbands,double);
	exp_w=ARRAY_2D(ncells,nbands,double);
	exp_temp=ARRAY_2D(ncells,nbands,double);
	
	mod_fmin=malloc(ncells*sizeof(double));
	mod_fmax=malloc(ncells*sizeof(double));
	J=malloc(ncells*sizeof(double));
	
	T_max=-1e99; //initialise the maximum temperature
	T_min=1e99; //initialise the maximum temperature
	
	
	for (i=0;i<ncells;i++)
	{
		mod_fmin[i]=1e99;
		mod_fmax[i]=-1e99;
		fscanf(contf, " %d", &idum1);
		fscanf(contf, " %d", &idum2);
		J[i]=0;
		for (j=0;j<nbands;j++)
		{
			fscanf(contf, "%le %le %d %le %le %le %le",&dum1,&dum2,&idum1,&dum3,&dum4,&dum5,&dum6);
			f1[i][j]=dum1;
			f2[i][j]=dum2;
			model[i][j]=idum1;
			pl_w[i][j]=dum3;
			pl_alpha[i][j]=dum4;
			exp_w[i][j]=dum5;
			exp_temp[i][j]=dum6;
			if (model[i][j]!=SPEC_MOD_FAIL && f1[i][j]<mod_fmin[i]) mod_fmin[i]=f1[i][j];
			if (model[i][j]!=SPEC_MOD_FAIL && f2[i][j]>mod_fmax[i]) mod_fmax[i]=f2[i][j];
			if (i==0) printf ("%e %e %d %e %e %e %e\n",f1[i][j],f2[i][j],model[i][j],pl_w[i][j],pl_alpha[i][j],exp_w[i][j],exp_temp[i][j]);
			if (model[i][j]!=SPEC_MOD_FAIL)
			{
				J[i]+=int_jnu(i,j);
			}
//			if (model[i][j]==SPEC_MOD_EXP) //Lines used in an attempt to speed things up with an exp lookup
//			{
//				if (exp_temp[i][j] > T_max) T_max=exp_temp[i][j];
//				if (exp_temp[i][j] < T_min) T_min=exp_temp[i][j];				
//			}

			
		}
	}
	
		
	printf ("Read spectral data - currently using %f Mb\n",g_usedMemory/1e6);
	printf ("In cell 0, spectral models run from %e to %e J=%e\n",mod_fmin[ncells-1],mod_fmax[ncells-1],J[ncells-1]);
	printf ("T_Max=%e T_min=%e\n",T_max,T_min);


	return(0);
}

/***************************************************************************/
//
// Read_line_data reads in all the levels and transitions.
//
/***************************************************************************/

int read_line_data()
{
	FILE *line_list, *atom_models;
	int z,istate,ion_index;
	int i,j,ib;
	int lev_tot;
	int idum1,idum2,idum3;
	double dum1,dum10,dum11,dum2;
	int nlvl,nlevels,nlvl_max;
	int nlines;
	int count;
	
    line_fmin=1e99;
    line_fmax=0.0;
	
	printf ("Reading atomic model data\n");	  
	//These lines look thruogh the atomic models file to see how much space we need to assign for models
	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open atomic_models.txt\n");
		abort();
	}
	nlvl_max=-1;	  
	while (fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1)==4)
	{
		for (i=0;i<nions;i++)
		{
			ion_index=-1;
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlevels[i]=nlevels=idum3;
				ion_info_xi[i]=dum1;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlevels=idum3;	
		
		for (nlvl = 0; nlvl < nlevels; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2); //index - level_energy - multiplicity - unused
 		}
	}
	fclose(atom_models);
	
	
	lev_energy=ARRAY_2D(nions,nlvl_max,double);
	lev_weight=ARRAY_2D(nions,nlvl_max,double);
	lev_pops=ARRAY_2D(nions,nlvl_max,double); //We will not be filling this array here - but lets set it up while we know nht max number of levels
	
	
	//Now read in the data for real
	
	count=0;
	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open adata.txt.\n");
		abort();
	}
	while (fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1)==4)
	{
		ion_index=-1;
		for (i=0;i<nions;i++) //Loop over ions to find the matching ion
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_index=i;
			}
			if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
				nlevels=idum3;
			else
			{
				nlevels=ion_info_nlevels[ion_index];
			}
		}
		for (nlvl = 0; nlvl < nlevels; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2); //index - level_energy - multiplicity - unused
			if (ion_index !=-1)
			{
				lev_weight[ion_index][nlvl]=dum11;
				lev_energy[ion_index][nlvl]=dum10;	
				
			}
 		}
	}
		

	fclose(atom_models);	
	printf ("Read level data - currently using %f Mb\n",g_usedMemory*1e-6);
		
	printf ("Reading line data\n");
	
	
	
	//These lines look thruogh the line data file to see how much space we need to assign for models
	if ((line_list = fopen("transitiondata.txt", "r")) == NULL)  //Try to open a file with transition data 
	{
		printf("Cannot open transitiondata.txt\n");
		abort();
	}
	nline_tot=0;
	while (fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3)==3)
	{
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlines[i]=nlines=idum3;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
				nline_tot=nline_tot+nlines;				
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlines=idum3;	
		for (nlvl = 0; nlvl < nlines; nlvl++)  //Loop over all the levels for this ion - we dont do anything with this data here
 		{
			fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1); //index - level_energy - multiplicity - unused
 		}
	}
	fclose(line_list);
	
	printf ("There are a total of %i lines linked to ions in our model\n",nline_tot);
	
	//Now we know how many lines we will be reading, we can set up the arrays
	
	
	trans_ion=malloc(nline_tot*sizeof(int));
	trans_lower=malloc(nline_tot*sizeof(int));
	trans_upper=malloc(nline_tot*sizeof(int));
	trans_iband=malloc(nline_tot*sizeof(int));
	
	
	trans_freq=malloc(nline_tot*sizeof(double));
	trans_lfreq=malloc(nline_tot*sizeof(double));
	trans_A_ul=malloc(nline_tot*sizeof(double));
	trans_B_ul=malloc(nline_tot*sizeof(double));
	trans_B_lu=malloc(nline_tot*sizeof(double));
	
	g_usedMemory+=4*nline_tot*sizeof(int);
	g_usedMemory+=5*nline_tot*sizeof(double);
	
	
	if ((line_list = fopen("transitiondata.txt", "r")) == NULL)  //Try to open a file with transition data 
	{
		printf("Cannot open transitiondata.txt\n");
		abort();
	}
	nline_tot=0;
	while (fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3)==3)
	{
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlines[i]=nlines=idum3;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlines=idum3;
		for (nlvl = 0; nlvl < nlines; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1); //index - level_energy - multiplicity - unused
			if (ion_index!=-1) //We only need to add this on if there is a matching ion in our  ion array
			{
				trans_ion[nline_tot]=ion_index;
				if (lev_energy[ion_index][idum3-1] > lev_energy[ion_index][idum2-1])
				{
					trans_lower[nline_tot]=idum2-1; //We subtract 1 because the level array starts at zero
					trans_upper[nline_tot]=idum3-1;
				}
				else //For some reason, the energy levels are swapped round - still use transition but swap the levels
				{
					trans_lower[nline_tot]=idum3-1;
					trans_upper[nline_tot]=idum2-1;
				}
				trans_freq[nline_tot]=(lev_energy[ion_index][trans_upper[nline_tot]]- lev_energy[ion_index][trans_lower[nline_tot]])*EV/H; //Compute the frequency for this transition
				trans_lfreq[nline_tot]=log10(trans_freq[nline_tot]);
	    		trans_A_ul[nline_tot]=dum1; //A_ul
	    		trans_B_ul[nline_tot]=CLIGHTSQUAREDOVERTWOH/pow(trans_freq[nline_tot],3.0) * trans_A_ul[nline_tot]; //B_ul
	    		trans_B_lu[nline_tot]=lev_weight[ion_index][trans_upper[nline_tot]]*trans_B_ul[nline_tot]/
															lev_weight[ion_index][trans_lower[nline_tot]]; //B_luq
				if (trans_freq[nline_tot]<line_fmin) line_fmin=trans_freq[nline_tot];
				if (trans_freq[nline_tot]>line_fmax) line_fmax=trans_freq[nline_tot];
				for (ib=0;ib<nbands;ib++) //Discover which band this line lies in - saves having to search every time in future.			
				{ 	
					if (trans_freq[nline_tot]>band_limits[ib] && trans_freq[nline_tot]<band_limits[ib+1])	trans_iband[nline_tot]=ib;
				}
				nline_tot++;
			}			
 		}
	}
	fclose(line_list);
	printf ("Fmin=%e Fmax=%e\n",line_fmin,line_fmax);
	printf ("Emin=%e Emax=%e\n",line_fmin/(EV/H),line_fmax/(EV/H));
	
	printf ("Read transition data - currently using %f Mb\n",g_usedMemory*1.e-6);	
	
	return(0);
}



/* ********************************************************************* */
char **Array2D (int nx, int ny, size_t dsize)
/*! 
 * Allocate memory for a 2-D array of any basic data type. - stolen from pluto!!!!
 *
 * \param [in] nx    number of elements in the 2nd dimension
 * \param [in] ny    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char **) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1]
 *          
 *********************************************************************** */
{
  int i;
  char **m;
 
  m    = (char **)malloc ((size_t) nx*sizeof(char *));
  m[0] = (char *) malloc ((size_t) nx*ny*dsize);
 
  for (i = 1; i < nx; i++) m[i] = m[(i - 1)] + ny*dsize;
 
  g_usedMemory += nx*ny*dsize;

  return m;
}


double model_jnu(double freq, double lfreq, int icell,int iband)
{
	int i,j;
	double j_nu;
	
	
	
	if (freq>f1[icell][iband] && freq<f2[icell][iband]) //The band could be truncated in this cell
	{
		if (model[icell][iband]==SPEC_MOD_FAIL)
			j_nu=0.0; //No model in this band - should have been caught already but belt and braces
		else if (model[icell][iband]==SPEC_MOD_PL)
			j_nu=pow (10, pl_w[icell][iband]+pl_alpha[icell][iband]*lfreq);
		else if (model[icell][iband]==SPEC_MOD_EXP)
			j_nu=exp_w[icell][iband] * exp1 ((-1.0 * H * freq) / (KB * exp_temp[icell][iband]));		
		else 
		{
			printf ("Serious problem getting a j_nu for freq %e\n",freq);
			abort();
		}				
	}
	else
	{
		j_nu=0.0; //There is no model - so no flux
	}
	
	return(j_nu);
	
		
}


double int_jnu(int icell,int iband)
{
	double a,b,integral,exp_pre;
		
	if (model[icell][iband]==SPEC_MOD_FAIL)
		integral=0.0; //No model in this band - should have been caught already but belt and braces
	else if (model[icell][iband]==SPEC_MOD_PL)
	{
	    a = pow (10.0, (pl_w[icell][iband] + (log10(f2[icell][iband]) * (pl_alpha[icell][iband] + 1.0))));
	    b = pow (10.0, (pl_w[icell][iband] + (log10(f1[icell][iband]) * (pl_alpha[icell][iband] + 1.0))));
		integral=(a - b) / (pl_alpha[icell][iband] + 1.0);
	}
	else if (model[icell][iband]==SPEC_MOD_EXP)
	{
	    exp_pre = (-1.0 * H) / (KB * exp_temp[icell][iband]);
		a =1./exp_pre*exp_w[icell][iband] * exp (exp_pre*f2[icell][iband]);
	    b =1./exp_pre*exp_w[icell][iband] * exp (exp_pre*f1[icell][iband]);
		integral=(a - b);		
	}
	else 
	{
		printf ("Serious problem integrating j_nu for band %i cell %i\n",iband,icell);
		abort();
	}
	
	return(integral);
}

double exp1(double x) {
	double ans;
	int indx;
	int test;
	
	if (fabs(x)<0.01)  //At this level - the difference between exp and the approximation is much less than 1pc
	{
		ans=1. + x*(1. + x/2.);
	}
	else if (fabs(x)>500)
	{
		ans=0.0;
	}
	else
	{ 
//		indx=(int)(-1000.*x+0.5);
//		ans=exp_lookup[indx-10];
//		printf ("%e %e %e\n",ans,exp(x),ans-exp(x));
		ans=exp(x);
	}
		
  
  return ans;
}
/* Attempt to speed things up with lookup - but not much faster and errors mounted up
double exp_lookup_init()
{
	int npoints,i;
	double exp_min,exp_max,test;
	double x,dx,xmin;
	dx=0.1;
	xmin=1.0;
	
	npoints=500000;
	
	exp_lookup=malloc(npoints*sizeof(double));
	
	for (i=0;i<npoints;i++)
	{
		x=(xmin+i*dx)/100.;
		exp_lookup[i]=exp(-1.*x);
	}
	
	

	return(npoints);
}
*/

  
