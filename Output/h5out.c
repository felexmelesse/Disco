
#include "../paul.h"
#include <hdf5.h>

void createFile( char * fname ){
   hid_t h5file = H5Fcreate( fname , H5F_ACC_TRUNC , H5P_DEFAULT , H5P_DEFAULT );
   H5Fclose( h5file );
}

void createGroup( char * fname , char * gname ){
   hid_t h5file = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gcreate1( h5file , gname , 0 );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void createDataset( char * fname , char * gname , char * dname , int dim , hsize_t * fdims , hid_t type ){
   hid_t h5file  = H5Fopen( fname , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5group = H5Gopen1( h5file , gname );
   hid_t fspace  = H5Screate_simple(dim,fdims,NULL);
   hid_t h5dset  = H5Dcreate1( h5group , dname , type , fspace , H5P_DEFAULT );
   H5Sclose( fspace );
   H5Dclose( h5dset );
   H5Gclose( h5group );
   H5Fclose( h5file );
}

void writeSimple( char * file , char * group , char * dset , void * data , hid_t type ){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );

   H5Dwrite( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

void writePatch( char * file , char * group , char * dset , void * data , hid_t type , int dim , int * start , int * loc_size , int * glo_size){
   hid_t h5fil = H5Fopen( file , H5F_ACC_RDWR , H5P_DEFAULT );
   hid_t h5grp = H5Gopen1( h5fil , group );
   hid_t h5dst = H5Dopen1( h5grp , dset );
   
   hsize_t mdims[dim];
   hsize_t fdims[dim];

   hsize_t fstart[dim];
   hsize_t fstride[dim];
   hsize_t fcount[dim];
   hsize_t fblock[dim];
   
   int d;
   for( d=0 ; d<dim ; ++d ){
      mdims[d] = loc_size[d];
      fdims[d] = glo_size[d];

      fstart[d]  = start[d];
      fstride[d] = 1;
      fcount[d]  = loc_size[d];
      fblock[d]  = 1;
   }
   hid_t mspace = H5Screate_simple(dim,mdims,NULL);
   hid_t fspace = H5Screate_simple(dim,fdims,NULL);

   H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , fblock );

   H5Dwrite( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

   H5Sclose( mspace );
   H5Sclose( fspace );
   H5Dclose( h5dst );
   H5Gclose( h5grp );
   H5Fclose( h5fil );
}

int Cell2Doub( struct cell * c , double * Q , int mode ){
   if( mode==0 ) return(NUM_Q+NUM_FACES+1); else{
      int q;
      for( q=0 ; q<NUM_Q ; ++q ) Q[q] = c->prim[q];
      for( q=0 ; q<NUM_FACES ; ++q ) Q[NUM_Q+q] = c->Phi[q];
      Q[NUM_Q+NUM_FACES] = c->piph;
      return(0);
   }
}

void dumpVal(char *filename, char *group, char *dset, void *val, 
                hid_t type)
{
    hsize_t fdims1[1];
    fdims1[0] = 1;
    createDataset(filename, group, dset, 1, fdims1, type);
    writeSimple(filename, group, dset, val, type);
}

void writeOpts(struct domain *theDomain, char filename[])
{
    char buf[256];
    char *buf2[1] = {buf};

    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, H5T_VARIABLE);

    createGroup(filename, "Opts");

    strncpy(buf, GIT_VERSION, 256);
    buf[255] = '\0';
    dumpVal(filename, "/", "GIT_VERSION", buf2, strtype);

    strncpy(buf, INITIAL, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "INITIAL", buf2, strtype);
    
    strncpy(buf, HYDRO, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "HYDRO", buf2, strtype);
    
    strncpy(buf, BOUNDARY, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "BOUNDARY", buf2, strtype);
    
    strncpy(buf, OUTPUT, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "OUTPUT", buf2, strtype);
    
    strncpy(buf, RESTART, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "RESTART", buf2, strtype);
    
    strncpy(buf, PLANETS, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "PLANETS", buf2, strtype);
    
    strncpy(buf, HLLD, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "HLLD", buf2, strtype);
    
    strncpy(buf, ANALYSIS, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "ANALYSIS", buf2, strtype);
    
    strncpy(buf, METRIC, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "METRIC", buf2, strtype);
    
    strncpy(buf, FRAME, 256);
    buf[255] = '\0';
    dumpVal(filename, "Opts", "FRAME", buf2, strtype);

    int val;
    
    val = NUM_C;
    dumpVal(filename, "Opts", "NUM_C", &val, H5T_NATIVE_INT);
    val = NUM_N;
    dumpVal(filename, "Opts", "NUM_N", &val, H5T_NATIVE_INT);
    val = CT_MODE;
    dumpVal(filename, "Opts", "CT_MODE", &val, H5T_NATIVE_INT);
}

void writePars(struct domain *theDomain, char filename[])
{
    struct param_list *pars = &(theDomain->theParList);

    createGroup(filename, "Pars");

    dumpVal(filename, "Pars", "Num_R",  &(pars->Num_R),  
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Num_Z",  &(pars->Num_Z),  
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "aspect", &(pars->aspect), 
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Num_Reports",  &(pars->NumRepts),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Num_Snapshots",  &(pars->NumSnaps),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Num_Checkpoints",  &(pars->NumChecks),
                    H5T_NATIVE_INT);

    dumpVal(filename, "Pars", "T_Start", &(pars->t_min),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "T_End", &(pars->t_max),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "R_Min", &(pars->rmin),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "R_Max", &(pars->rmax),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Z_Min", &(pars->zmin),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Z_Max", &(pars->zmax),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Phi_Max", &(pars->phimax),
                    H5T_NATIVE_DOUBLE);
    
    dumpVal(filename, "Pars", "Use_Logtime", &(pars->Out_LogTime),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Log_Zoning", &(pars->LogZoning),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Z_Periodic", &(pars->Z_Periodic),
                    H5T_NATIVE_INT);

    dumpVal(filename, "Pars", "Log_Radius", &(pars->LogRadius),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Max_Aspect_Short", &(pars->MaxShort),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Max_Aspect_Long", &(pars->MaxLong),
                    H5T_NATIVE_DOUBLE);

    dumpVal(filename, "Pars", "CFL", &(pars->CFL),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "PLM", &(pars->PLM),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Max_DT", &(pars->maxDT),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Adiabatic_Index", &(pars->Adiabatic_Index),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Isothermal", &(pars->isothermal_flag),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Cs2_Profile", &(pars->Cs2_Profile),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Cs2_Par", &(pars->Cs2_Par),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Density_Floor", &(pars->Density_Floor),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Pressure_Floor", &(pars->Pressure_Floor),
                    H5T_NATIVE_DOUBLE);

    dumpVal(filename, "Pars", "Mesh_Motion", &(pars->Mesh_Motion),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Riemann_Solver", &(pars->Riemann_Solver),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Initial_Regrid", &(pars->Initial_Regrid),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Restart", &(pars->restart_flag),
                    H5T_NATIVE_INT);

    dumpVal(filename, "Pars", "Exact_Mesh_Omega", &(pars->Exact_Mesh_Omega),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Exact_Mesh_Omega_Par", 
                    &(pars->Exact_Mesh_Omega_Par), H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Energy_Omega", &(pars->Energy_Omega),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Energy_Omega_Par", 
                    &(pars->Energy_Omega_Par), H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "RotFrame", &(pars->RotFrame),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "RotOmega", &(pars->RotOmega),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "RotD", &(pars->RotD),
                    H5T_NATIVE_DOUBLE);
    
    dumpVal(filename, "Pars", "Use_Viscosity", &(pars->visc_flag),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Viscosity", &(pars->viscosity),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Use_As_Alpha", &(pars->alpha_flag),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Include_Atmos", &(pars->include_atmos),
                    H5T_NATIVE_INT);
    
    dumpVal(filename, "Pars", "Mach_Number", &(pars->Disk_Mach),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Mass_Ratio", &(pars->Mass_Ratio),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Eccentricity", &(pars->Eccentricity),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Drift_Rate", &(pars->Drift_Rate),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Drift_Exp", &(pars->Drift_Exp),
                    H5T_NATIVE_DOUBLE);
    
    dumpVal(filename, "Pars", "Constrained_Transport", &(pars->CT),
                    H5T_NATIVE_INT);
    
    dumpVal(filename, "Pars", "Metric_Par0", &(pars->metricPar0),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Metric_Par1", &(pars->metricPar1),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Metric_Par2", &(pars->metricPar2),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Metric_Par3", &(pars->metricPar3),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Metric_Par4", &(pars->metricPar4),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Init_Par0", &(pars->initPar0),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Init_Par1", &(pars->initPar1),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Init_Par2", &(pars->initPar2),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Init_Par3", &(pars->initPar3),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Init_Par4", &(pars->initPar4),
                    H5T_NATIVE_DOUBLE);
    
    dumpVal(filename, "Pars", "Noise_Type", &(pars->noiseType),
                    H5T_NATIVE_INT);
    dumpVal(filename, "Pars", "Noise_Abs", &(pars->noiseAbs),
                    H5T_NATIVE_DOUBLE);
    dumpVal(filename, "Pars", "Noise_Rel", &(pars->noiseRel),
                    H5T_NATIVE_DOUBLE);
}

double get_dV( double * , double * );
void prim2cons( double * , double * , double *, double );
void cons2prim( double * , double * , double *, double );

void zero_diagnostics( struct domain * );
void avg_diagnostics( struct domain * );

void output( struct domain * theDomain , char * filestart ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   int Ng = theDomain->Ng;
   int Nr_Tot = theDomain->theParList.Num_R;
   int Nz_Tot = theDomain->theParList.Num_Z;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Z_Periodic = theDomain->theParList.Z_Periodic;

   int NpDat = 6;
   int Ntools = theDomain->num_tools;
   avg_diagnostics( theDomain );

   char filename[256];
   sprintf(filename,"%s.h5",filestart);
   int jmin = 0;
   if( dim_rank[0] != 0 ) jmin = Ng;
   int jmax = Nr;
   if( dim_rank[0] != dim_size[0]-1 ) jmax = Nr-Ng;
   int kmin = 0;
   //if( dim_rank[1] != 0 || Z_Periodic ) kmin = Ng;
   if( dim_rank[1] != 0 ) kmin = Ng;
   int kmax = Nz;
   //if( dim_rank[1] != dim_size[1]-1 || Z_Periodic ) kmax = Nz-Ng;
   if( dim_rank[1] != dim_size[1]-1 ) kmax = Nz-Ng;

   if(Z_Periodic)
   {
       Nz_Tot += 2*Ng;
   }

   int Ntot = 0;
   int j,k;
   for( j=jmin ; j<jmax ; ++j ){
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         Ntot += Np[jk];
      }
   }
   int myNtot = Ntot;
   MPI_Allreduce( MPI_IN_PLACE , &Ntot  , 1 , MPI_INT , MPI_SUM , theDomain->theComm );

   int Ndoub = Cell2Doub(NULL,NULL,0);

   hsize_t fdims1[1];
   hsize_t fdims2[2];
   if( rank==0 ){
      printf("Writing Checkpoint...\n");
      
      createFile(filename);
      createGroup(filename,"Grid");

      fdims1[0] = 1;
      createDataset(filename,"Grid","T",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Nr_Tot+1;
      createDataset(filename,"Grid","r_jph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims1[0] = Nz_Tot+1;
      createDataset(filename,"Grid","z_kph",1,fdims1,H5T_NATIVE_DOUBLE);
      fdims2[0] = Nz_Tot;
      fdims2[1] = Nr_Tot;
      createDataset(filename,"Grid","Index",2,fdims2,H5T_NATIVE_INT);
      createDataset(filename,"Grid","Np",2,fdims2,H5T_NATIVE_INT);
      createDataset(filename,"Grid","Id_phi0",2,fdims2,H5T_NATIVE_INT);

      createGroup(filename,"Data");

      fdims2[0] = Ntot;
      fdims2[1] = Ndoub;
      createDataset(filename,"Data","Cells",2,fdims2,H5T_NATIVE_DOUBLE);
      fdims2[0] = Npl;
      fdims2[1] = NpDat;
      createDataset(filename,"Data","Planets",2,fdims2,H5T_NATIVE_DOUBLE);
      hsize_t fdims3[3] = {Nz_Tot, Nr_Tot, Ntools};
      createDataset(filename,"Data","Diagnostics",3,fdims3,H5T_NATIVE_DOUBLE);
   }
   MPI_Barrier( theDomain->theComm );
   if( rank==0 ){
      writeSimple(filename,"Grid","T",&(theDomain->t),H5T_NATIVE_DOUBLE);
      writePars(theDomain, filename);
      writeOpts(theDomain, filename);
      double PlanetData[Npl*NpDat];
      int p;
      for( p=0 ; p<Npl ; ++p ){
         struct planet * pl = theDomain->thePlanets+p;
         PlanetData[NpDat*p + 0] = pl->M;
         PlanetData[NpDat*p + 1] = pl->vr;
         PlanetData[NpDat*p + 2] = pl->omega;
         PlanetData[NpDat*p + 3] = pl->r;
         PlanetData[NpDat*p + 4] = pl->phi;
         PlanetData[NpDat*p + 5] = pl->eps;
      }
      writeSimple(filename,"Data","Planets",PlanetData,H5T_NATIVE_DOUBLE);
   }

   int jSize = jmax-jmin;
   int kSize = kmax-kmin;
   int nrk;
   int j0 = -1;
   int k0 = -1;
   int jSum = 0;
   for( nrk=0 ; nrk < dim_size[0] ; ++nrk ){
      if( nrk == dim_rank[0] ){
         j0 = jSum;
         if( dim_rank[1] == 0 ) jSum += jSize;
      }
      MPI_Allreduce( MPI_IN_PLACE , &jSum , 1 , MPI_INT , MPI_MAX , theDomain->theComm );
   }
   int kSum = 0;
   for( nrk=0 ; nrk < dim_size[1] ; ++nrk ){
      if( nrk == dim_rank[1] ){
         k0 = kSum;
         if( dim_rank[0] == 0 ) kSum += kSize;
      }
      MPI_Allreduce( MPI_IN_PLACE , &kSum , 1 , MPI_INT , MPI_MAX , theDomain->theComm );
   }


   if( Nr_Tot == 1 ){ j0 = 0; jSum = 1; }
   if( Nz_Tot == 1 ){ k0 = 0; kSum = 1; }


   int * Index   = (int *) malloc( jSize*kSize*sizeof(int) );
   int * Size    = (int *) malloc( jSize*kSize*sizeof(int) );
   int * Id_phi0 = (int *) malloc( jSize*kSize*sizeof(int) );
   double * diagRZwrite = (double *) malloc( jSize*kSize*Ntools*sizeof(double) );
   double * Qwrite = (double *) malloc( myNtot*Ndoub*sizeof(double) );

   double *Qrz = theDomain->theTools.Qrz;


   int index = 0;
   int q;
   for( k=kmin ; k<kmax ; ++k ){
      for( j=jmin ; j<jmax ; ++j ){
         int jk = (k-kmin)*jSize + (j-jmin);
         Index[jk] = index;
         Size[jk] = Np[j+Nr*k];
         for(q=0; q<Ntools; q++)
            diagRZwrite[Ntools*jk+q] = Qrz[Ntools*(j+Nr*k)+q];
 
         double phi0 = M_PI;
         int Id = 0;
         int i;
         for( i=0 ; i<Np[j+Nr*k] ; ++i ){
            struct cell * c = &(theCells[j+Nr*k][i]);
            Cell2Doub( c , Qwrite+index*Ndoub , 1 );
            double phi = c->piph-.5*c->dphi;
            if( cos(phi0) < cos(phi) ){ phi0 = phi; Id = index; }
            ++index;
         }
         Id_phi0[jk] = Id;
      }
   }

   int runningTot = 0;
   for( nrk=0 ; nrk < size ; ++nrk ){
      int thisTot = myNtot;
      MPI_Bcast( &thisTot , 1 , MPI_INT , nrk , theDomain->theComm );
      if( rank > nrk ) runningTot += thisTot;
   }

   int jk;
   for( jk=0 ; jk<jSize*kSize ; ++jk ){
      Index[jk] += runningTot;
      Id_phi0[jk] += runningTot;
   }

   for( nrk=0 ; nrk < size ; ++nrk ){
      if( rank==nrk ){      
         //Write Cell Data
         int start2[2]    = {runningTot,0};
         int loc_size2[2] = {myNtot,Ndoub};
         int glo_size2[2] = {Ntot,Ndoub};
         writePatch( filename , "Data" , "Cells" , Qwrite , H5T_NATIVE_DOUBLE , 2 , start2 , loc_size2 , glo_size2 );
         //Write Indices and Sizes for each radial track
         start2[0] = k0;
         start2[1] = j0;
         loc_size2[0] = kSize;
         loc_size2[1] = jSize;
         glo_size2[0] = Nz_Tot;
         glo_size2[1] = Nr_Tot;
         writePatch( filename , "Grid" , "Index"   , Index   , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         writePatch( filename , "Grid" , "Np"      , Size    , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         writePatch( filename , "Grid" , "Id_phi0" , Id_phi0 , H5T_NATIVE_INT , 2 , start2 , loc_size2 , glo_size2 );
         int start3[3] = {k0, j0, 0};
         int loc_size3[3] = {kSize, jSize, Ntools};
         int glo_size3[3] = {Nz_Tot, Nr_Tot, Ntools};
         writePatch( filename , "Data" , "Diagnostics" , diagRZwrite , H5T_NATIVE_DOUBLE , 3 , start3 , loc_size3 , glo_size3 );

         //Write 1D Radial Data
         if( dim_rank[1] == 0 ){
            int offset = Ng;
            if( dim_rank[0] == 0 ) offset = 0;
            int start1[1]    = {j0};
            int loc_size1[1] = {jSize};
            if( dim_rank[0] == dim_size[0]-1 ) loc_size1[0]++;
            int glo_size1[1] = {Nr_Tot+1};
            writePatch( filename , "Grid" , "r_jph" , r_jph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
         }
         //Write 1D Vertical Data
         if( dim_rank[0] == 0 ){
            int offset = Ng;
            if( dim_rank[1] == 0 ) offset = 0;
            //if( Z_Periodic && dim_rank[1] == 0 ) offset += Ng;
            int start1[1]    = {k0};
            int loc_size1[1] = {kSize};
            if( dim_rank[1] == dim_size[1]-1 ) loc_size1[0]++;
            int glo_size1[1] = {Nz_Tot+1};
            writePatch( filename , "Grid" , "z_kph" , z_kph-1+offset , H5T_NATIVE_DOUBLE , 1 , start1 , loc_size1 , glo_size1 );
         }
      }
      MPI_Barrier( theDomain->theComm );
   }
   zero_diagnostics( theDomain );

   free(Index);
   free(Size);
   free(Id_phi0);
   free(Qwrite);
   free(diagRZwrite);
   MPI_Barrier(theDomain->theComm);
}


