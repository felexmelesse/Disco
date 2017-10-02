#include "paul.h"

static int Nmc = 0.0;

double get_dV( double *, double * );

void setTracerParams( struct domain * theDomain){

   int initType = theDomain->theParList.tr_init_type;
   int num_tracers = theDomain->theParList.num_tracers;

   if( !initType ){
      int Ncells = theDomain->Ncells;
      theDomain->Ntr = Ncells;
   }else if( initType==1 ){
      if( theDomain->rank==0 )
         printf("Total Tracers: %d \n", num_tracers );

      double rmin = theDomain->theParList.rmin;
      double rmax = theDomain->theParList.rmax;
      double zmin = theDomain->theParList.zmin;
      double zmax = theDomain->theParList.zmax;
      double phi_max = theDomain->theParList.phimax;
      double XM[3] = {rmin, 0.0, zmin};
      double XP[3] = {rmax, phi_max, zmax};
      double Vtot = get_dV( XP, XM );
      double ratio = num_tracers/Vtot;

      double r0   = theDomain->r0;
      double z0   = theDomain->z0;
      double delr = theDomain->delr;
      double delz = theDomain->delz;
      double xm[3] = {r0, 0.0, z0};
      double xp[3] = {r0+delr, phi_max, z0+delz};
      double dV = get_dV( xp, xm );

      theDomain->Ntr = dV*ratio;
   }else{
      int size = theDomain->size;
      theDomain->Ntr = num_tracers/size;
   }

   theDomain->Nmc = Nmc;
   //theDomain->tr_out = theDomain->theParList.tr_out_flag;
   MPI_Barrier( theDomain->theComm );
   printf("Rank %d: Ntr = %d \n", theDomain->rank, theDomain->Ntr);

}

double getRandIn( double xmin, double dx){

  return xmin + ( (double)rand()/(double)RAND_MAX )*dx;
}

void printTracerCoords( struct domain * );

void initTracers_Rand( struct domain *theDomain ){  //randomly init tracers in serial (whole domain)

   struct param_list theParamList = theDomain->theParList;
   double rmin = theParamList.rmin;
   double rmax = theParamList.rmax;
   double dr   = rmax - rmin;
   double zmin = theParamList.zmin;
   double zmax = theParamList.zmax;
   double dz   = zmax - zmin;
   double phimax = theParamList.phimax;

   srand(theDomain->rank);
   rand();

   struct tracerList *theList = theDomain->theTracers;
   struct tracer *tr = theList->head;
   while( tr!=NULL ){
     double r = getRandIn( rmin, dr );
     double z = getRandIn( zmin, dz );
     double phi = getRandIn( 0.0, phimax );
     tr->R = r; tr->Z = z; tr->Phi = phi;
     tr->Type   = 1;
     tr->rmFlag = 0;
     tr = tr->next;
   }
   //printTracerCoords( theDomain );
}


int getN0( int , int , int );
double get_moment_arm( double *, double * );
void syncTracerIDs( struct domain * );

void initializeTracers( struct domain *theDomain ){

  int initType = theDomain->theParList.tr_init_type;
  int id = 0.0;
  int rank = theDomain->rank;
  struct tracer *tr = theDomain->theTracers->head;

  if( !initType ){
      if( !rank )
         printf("Initializing tracers per cell\n");

      struct cell **theCells = theDomain->theCells;
      int  Nr = theDomain->Nr;
      int  Nz = theDomain->Nz;
      int *Np = theDomain->Np;
      double *r_jph = theDomain->r_jph;
      double *z_kph = theDomain->z_kph;

      int i,j,k;
      for( j=0; j<Nr; ++j ){
         for( k=0; k<Nz; ++k ){
            int jk = j+Nr*k;
            double rm = r_jph[j-1];
            double rp = r_jph[j];
            double zm = z_kph[k-1];
            double zp = z_kph[k];
            for( i=0; i<Np[jk]; ++i ){
               if( tr==NULL ){
                  printf("Tracer list size does not match cell number\n");
                  break;
               }
               struct cell *c = &(theCells[jk][i]);
               double phip = c->piph;
               double phim = phip - c->dphi;
               double xp[3] = {rp,phip,zp};
               double xm[3] = {rm,phim,zm};
               double r = get_moment_arm(xp,xm);
               tr->R   = r;
               tr->Phi = 0.5*(phim+phip);
               tr->Z   = 0.5*(zm+zp);
               tr->ID  = id;
               tr->Type   = 1;
               tr->Charac = 0.0;
               tr->rmFlag = 0;
               id++;
               tr = tr->next;
            }
         }
      }
  }else{
      if( !rank )
         printf("Initializing tracers by number\n");

      double r0   = theDomain->r0;
      double z0   = theDomain->z0;
      double delr = theDomain->delr;
      double delz = theDomain->delz;
      double phi_max = theDomain->theParList.phimax;

      srand(rank);
      rand();

      double r, z, phi;
      while( tr != NULL){
         r = getRandIn(r0, delr);
         z = getRandIn(z0, delz);
         phi = getRandIn(0, phi_max);
         tr->R = r; tr->Z = z; tr->Phi = phi;
         tr->ID     = id;
         tr->Type   = 1; // rank+1;
         tr->rmFlag = 0;
         id++;
         tr = tr->next;
     }
  }
  //printTracerCoords( theDomain );
  syncTracerIDs( theDomain );
}


int check_phi(double phi, double phip, double dphi, double phi_max){

   double phic = phi - phip;
   while( phic >  phi_max/2.0 ){ phic -= phi_max; }
   while( phic < -phi_max/2.0 ){ phic += phi_max; }
   if( -dphi<phic && phic<0 ){ return 1; }
   else{ return 0; }

}

int check_in_cell(struct tracer *tr, double *xp, double *xm, double phi_max){

   double r   = tr->R;
   double phi = tr->Phi;
   double z   = tr->Z;

   if( xm[0] < r && r < xp[0] ){
	if( xm[2] < z && z < xp[2]){
	     if( check_phi(phi, xp[1], xp[1]-xm[1], phi_max) ){
   		return 1;
	     }
	}
   }
   return 0;
}

void test_cell_vel( struct tracer *tr, struct cell *c ){

   int check1=0, check2=0, check3=0;
   if( isnan(c->prim[URR]) ){
	   tr->Vr 	  = 0;
	   tr->Omega = 0;
   	tr->Vz    = 0;
	   //tr->Type  = 3;
   	check1 = 1;
   }
   if( isnan(c->prim[UPP]) ){
	   tr->Vr    = 0;
      tr->Omega = 0;
      tr->Vz    = 0;
      //tr->Type  = 4;
   	check2 = 1;
   }
   if( isnan(c->prim[UZZ]) ){
	   tr->Vr    = 0;
      tr->Omega = 0;
      tr->Vz    = 0;
      //tr->Type  = 5;
	   check3 = 1;
   }

  if( check1==1 && check2==1 && check3==1){ tr->Type = 6; }

}


void get_local_vel(struct tracer *tr, struct cell *c){

   double vr, om, vz;
   int type = 1;

   //For testing
   double Vv  = sqrt(20);
   double phi = tr->Phi;

  if( c != NULL ){
      vr = c->prim[URR];
      om = c->prim[UPP];
      vz = c->prim[UZZ];
   } else{
      vr = 0.0;
      om = 0.0;
      vz = 0.0;
      type = 2;
      tr->Type = type;
   }
   tr->Vr    = vr;
   tr->Omega = om;
   tr->Vz    = vz;
   //tr->Type  = type;
   //if( c != NULL ){
   //  test_cell_vel( tr, c );
   //}
}

void set_tracer_vel( struct tracer *tr ){
   //For testing
   double Vv  = sqrt(50);
   double phi = tr->Phi;

   tr->Vr   =  Vv*cos(phi);
   tr->Omega= -Vv*sin(phi);
   tr->Vz   =  Vv*Vv;

}

int binSearch(double target, double *arr, int size){

  int lo = 0;
  int hi = size;
  int mid = lo + (hi-lo)/2;
  while( lo < hi ){
    if( arr[mid] < target )
      lo = mid+1;
    else
      hi = mid;
    mid = lo + (hi-lo)/2;
  }
  return mid;
}

struct cell * get_tracer_cell(struct domain *theDomain, struct tracer *tr){

   struct cell **theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   double phi_max = theDomain->phi_max;
   double *r_jph = theDomain->r_jph;
   double *z_kph = theDomain->z_kph;

   int i;
   int j  = binSearch( tr->R, r_jph, Nr );
   int k  = binSearch( tr->Z, z_kph, Nz );
   int jk = j + Nr*k;
   for( i=0; i<Np[jk]; ++i ){
     struct cell *c = &(theCells[jk][i]);
     if( check_phi( tr->Phi, c->piph, c->dphi, phi_max) )
        return c;
   }
/*
   int i,j,k;
   for( j=0; j<Nr; ++j){
	for( k=0; k<Nz; ++k){
	   int jk = j+Nr*k;
	   double rm = r_jph[j-1];
	   double rp = r_jph[j];
	   double zm = z_kph[k-1];
	   double zp = z_kph[k];
	   for( i=0; i<Np[jk]; ++i){
	   	struct cell *c = &(theCells[jk][i]);
	   	double phip = c->piph;
		   double phim = phip - c->dphi;
		   double xp[3] = {rp, phip, zp};
		   double xm[3] = {rm, phim, zm};
		   if( check_in_cell(tr, xp, xm, phi_max) ){
         printf("(i,j,k): (%d, %d, %d)\n", i,j,k);
         printf("(i,j,k)_1: (%d, %d, %d)\n", i1, j1, k1);
			   return c;
		}
	   }
	}
   }
*/
   return NULL;

}


void moveTracers(struct domain *theDomain, struct tracer *tr, double dt){

   double rmin = theDomain->theParList.rmin;
   double rmax = theDomain->theParList.rmax;
   double zmin = theDomain->theParList.zmin;
   double zmax = theDomain->theParList.zmax;

   double r = tr->R;
   double phi = tr->Phi;
   double z = tr->Z;

   r += tr->Vr*dt;
   if( r > rmax ){
      r = 0.0/0.0;
      //tr->Type = 11;
   }
   if( r < rmin ) {
      r = 0.0/0.0;
      //tr->Type = 11;
   }

   z += tr->Vz*dt;
   if( z > zmax ){
      z = 0.0/0.0;
      //tr->Type = 21;
   }
   if( z < zmin ){
      z = 0.0/0.0;
      //tr->Type = 21;
   }

   tr->R = r;
   tr->Z = z;

   phi += tr->Omega*dt;
   tr->Phi = phi;
}


void tracer_RK_copy( struct tracer * tr ){
   tr->RK_r    = tr->R;
   tr->RK_phi  = tr->Phi;
   tr->RK_z	   = tr->Z;
   tr->RK_vr	= tr->Vr;
   tr->RK_omega= tr->Omega;
   tr->RK_vz	= tr->Vz;
}

void tracer_RK_adjust( struct tracer * tr , double RK ){
   tr->R     = (1.-RK)*tr->R     + RK*tr->RK_r;
   tr->Phi   = (1.-RK)*tr->Phi   + RK*tr->RK_phi;
   tr->Z     = (1.-RK)*tr->Z     + RK*tr->RK_z;
   tr->Vr    = (1.-RK)*tr->Vr    + RK*tr->RK_vr;
   tr->Omega = (1.-RK)*tr->Omega + RK*tr->RK_omega;
   tr->Vz    = (1.-RK)*tr->Vz    + RK*tr->RK_vz;
}


void updateTracers(struct domain *theDomain, double dt){

   //printf("Trying to update tracers!");
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
   	struct cell   *c  = get_tracer_cell( theDomain , tr );
	   get_local_vel( tr , c );
      //set_tracer_vel( tr );
      //if( theDomain->rank==1 ) printf("%f\n", (tr->Vr)*(tr->Vr)+(tr->Omega)*(tr->Omega) );
      moveTracers( theDomain , tr , dt );
      tr = tr->next;
   }

}
