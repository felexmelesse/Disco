#include "paul.h"

static int      Nmc = 0;
static double   tr_rmax = 0.0;
static double   tr_rmin = 0.0;
static double   BIG_NUM = 10000;
static double   gamma_law = 0.0;
//static double r_mini = 1.0;
//static double r_gap  = 2.0;

void set_tracer_bounds( double *, double *);
double SS_area( double, double );
double get_dV( double *, double * );

void setTracerParams( struct domain * theDomain){

   //r_sink = theDomain->theParList.r_sink;
   gamma_law = theDomain->theParList.Adiabatic_Index;
   int initType = theDomain->theParList.tr_init_type;
   int num_tracers = theDomain->theParList.num_tracers;
   tr_rmax = theDomain->theParList.tr_rmax;
   tr_rmin = theDomain->theParList.tr_rmin;
   if( tr_rmin > tr_rmax )
       printf("WARNING: tracer r_min > r_max \n");

   if( !initType ){
     //number of tracers = number of cells--doesn't have tr_max yet
      int Ncells = theDomain->Ncells;
      theDomain->Ntr = Ncells;
   }else if( initType==1 || initType==3 ){
      if( theDomain->rank==0 )
         printf("Total Tracers: %d \n", num_tracers );

      double rmin = theDomain->theParList.rmin;
      double rmax = theDomain->theParList.rmax;
      if( tr_rmax > 0.0 && rmax > tr_rmax )
          rmax = tr_rmax;
      if( tr_rmin > 0.0 && tr_rmin > rmin )
          rmin = tr_rmin;
      double zmin = theDomain->theParList.zmin;
      double zmax = theDomain->theParList.zmax;
      double phi_max = theDomain->theParList.phimax;
      double XM[3] = {rmin, 0.0, zmin};
      double XP[3] = {rmax, phi_max, zmax};
      
      double r0   = theDomain->r0;
      double z0   = theDomain->z0;
      double delr = theDomain->delr;
      double delz = theDomain->delz;
      double r1   = r0 + delr;

       
      //set_tracer_bounds( &r0, &r1 );
      //TODO: Move these too some function that sets r1 and r0 only for 
      //      tracer initialization
      //if( tr_rmax > 0.0 ){
      
      if( tr_rmax ){
        if( tr_rmax > r0 ){
            if( r1 > tr_rmax )
                r1 = tr_rmax;
        }
        else{
            r1 = r0; //No tracers in this process
        }
      }
      if( tr_rmin ){
        if( tr_rmin < r1 ){
            if( tr_rmin > r0 )
                r0 = tr_rmin;
        }
        else{
            r1 = r0; //No tracers in this process
        }
      }

      double xm[3] = {r0, 0.0, z0};
      double xp[3] = {r1, phi_max, z0+delz};
      if( initType==1 ){
          double Vtot = get_dV( XP, XM );
          double ratio = num_tracers/Vtot;
          double dV = get_dV( xp, xm );
          theDomain->Ntr = dV*ratio;
      }
      else{
          double ratio = num_tracers/SS_area( rmin, rmax );
          theDomain->Ntr = ratio*SS_area( r0, r1 );
      }
   }else{
      //equal number of tracers by process--doesn't have tr_rmax yet
      int size = theDomain->size;
      theDomain->Ntr = num_tracers/size;
   }

   theDomain->Nmc = Nmc;
   //theDomain->tr_out = theDomain->theParList.tr_out_flag;
   MPI_Barrier( theDomain->theComm );
   //printf("Rank %d: Ntr = %d \n", theDomain->rank, theDomain->Ntr);

}

double SS_area( double r1, double r2 ){
    //Area under the curve r^(-1/2) from r1 to r2
    return 2*( sqrt(r2) - sqrt(r1) );
}

void set_tracer_bounds( double* r0, double* r1 ){

    if( tr_rmax ){
        if( tr_rmax > (*r0) ){
            if( (*r1) > tr_rmax )
                (*r1) = tr_rmax;
        }
        else{
            (*r1) = (*r0); //No tracers in this process
        }
    }
    if( tr_rmin ){
        if( tr_rmin < (*r1) ){
            if( tr_rmin > (*r0) )
                (*r0) = tr_rmin;
        }
        else{
            (*r1) = (*r0); //No tracers in this process
        }
    }
    return;
}

double getRandIn( double xmin, double dx){

  return xmin + ( (double)rand()/(double)RAND_MAX )*dx;
}

double getRand_SS( double xmin, double norm_fac ){

    double u  = (double)rand()/(double)RAND_MAX;
    double rs = 0.5*u*norm_fac + sqrt(xmin);
    return rs*rs;
    //return xmin + ( (double)rand()/(double)RAND_MAX )*dx;
}

/*
void tr_gap_split( struct domain *theDomain ){

   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      double r = tr->R;
      if( r < r_gap )
         tr->Type = 50;
      tr = tr->next;
   }
}
*/

void printTracerCoords( struct domain * );
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
               tr->myCell = c;
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
      
      double r1 = r0 + delr;
      if( tr_rmax > 0.0 && r1 > tr_rmax ){
          r1 = tr_rmax;
          delr = r1 - r0;
          //delr will be negative for processes after tr_rmax
          //shouldn't matter bc should get 0 tracers
      }
      if( tr_rmin > 0.0 && r0 < tr_rmin ){
          r0 = tr_rmin;
          delr = r1 - r0;
      }
      srand(rank);
      rand();

      //double r_cav = 2.0;
      double r, z, phi;
      while( tr != NULL){
         if( initType == 3 ){
            double norm_fac = SS_area( r0, r1 );
            r = getRand_SS(r0, norm_fac);
         }
         else
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
  //tr_gap_split( theDomain );
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


void get_local_vel( struct tracer *tr ){

    struct cell *c = tr->myCell;
    double vr, om, vz;
    int type = 1;

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

void get_tracer_cell(struct domain *theDomain, struct tracer *tr){

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
     if( check_phi( tr->Phi, c->piph, c->dphi, phi_max) ){
        tr->myCell = c;
        //set dPdphi
        //tr->dPdp = c->gradp[PPP];
        //if( tr->dPdp >  BIG_NUM )
        //    tr->dPdp = 0.0/0.0;
        return;
     }
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
   tr->myCell = NULL;
   return;

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
   if( r > rmax )
      r = 0.0/0.0;
   if( r < rmin )
      r = 0.0/0.0;

   z += tr->Vz*dt;
   if( z > zmax )
      z = 0.0/0.0;
   if( z < zmin )
      z = 0.0/0.0;

   tr->R = r;
   tr->Z = z;

   phi += tr->Omega*dt;
   tr->Phi = phi;

/*
   if(r < r_mini && tr->Type==1 )
      tr->Type = 99;
*/
}

double get_dvdt( struct tracer *tr ){
    
    double dvr = tr->Vr - tr->v_prev[0];
    double dvp = (tr->Omega)*(tr->R) - tr->v_prev[0];
    double dvz = tr->Vz - tr->v_prev[2];

    return sqrt( dvr*dvr + dvp*dvp + dvz*dvz );
}   

//Specific to Binary runs -> at present doing this for expediency
//double get_planets_phi( double , double );
double get_potential( double, double );
double get_cs2( double, double );

void updatePhysQuants( struct tracer * tr, double gamma ){
    //Need to update (and treatment of tracer charactersitics
    //overall) to be more general
    struct cell *c = tr->myCell;
    double rho_0 = tr->prim[RHO];
    int q;
    if( c!=NULL ){
        //give prim vars to tracer
        for( q=0; q<NUM_Q; q++ ){
            double x = c->prim[q];
            if( x > BIG_NUM ) 
                x = -1.0;
            tr->prim[q] = x;
        }    
    }
    else{
        for( q=0; q<NUM_Q; q++ )
            tr->prim[q] = -1.0;
    }
        
//================ Calculate jacobi constant ================================
    double r  = tr->R;
    double vp = r*(tr->Omega - 1.0 ); //Make sure velocity is in corot frame
    double v2 = (tr->Vr)*(tr->Vr) + vp*vp + (tr->Vz)*(tr->Vz);
//    double bin_pot = get_planets_phi( r , tr->Phi );
    double bin_pot = get_potential( r , tr->Phi );
    double C_j = r*r - v2 + 2*bin_pot;
    if( fabs(C_j) > BIG_NUM )
        C_j = -1;
    tr->C_j = C_j;
//===========================================================================
//
//================ Calculate Change in C_j ==================================
//                 based on D.J. D'Orazio et al. (2016)
    double rho = tr->prim[RHO];
    double cs2 = get_cs2( r, tr->Phi );
    double ratio = rho/rho_0;
    if( ratio < 0.0 ) //maybe change how this is handled
        ratio = 1.0;
    double dCj = cs2/gamma_law*log( ratio ); // = /int dP/rho = c_s^2/gamma *ln(rho/rho_old)
    if( fabs(dCj) > BIG_NUM )
        dCj = 0.0;
    tr->dPdp = dCj;
//===========================================================================


/*    
    double r = tr->R;
    double Lp = -1.0;
    double Jj = -1.0;
    double Ss = -1.0;
    double rho = 1.0;
    if( c!=NULL ){
        rho = c->prim[RHO];
        //double rho = c->prim[RHO]
        double Pp  = c->prim[PPP];
        double om  = c->prim[UPP];
        double r2  = r*r;
        Lp = om*r2;   // *rho;
        Ss = log( Pp/pow(rho,gamma) );
        Jj = Lp - pow(r,-1.5)*r2; 
    }
    if( Lp > BIG_NUM ) Lp = -1.0;
    if( Ss > BIG_NUM ) Ss = -1.0;
    if( Jj > BIG_NUM ) Jj = -1.0;
    
    tr->Lp = Lp*rho;
    tr->Jj = Jj;         //specific angmom minus keplerian background
    tr->Ss = Ss;
*/
/*
    double fgrav = -1;
    double gravfrac = -1;
    if( c!=NULL ){
        //give prim vars to tracer
        //int q = 0;
        for( q=0; q<NUM_Q; q++ ){
            double x = c->prim[q];
            if( x > BIG_NUM ) 
                x = -1.0;
            tr->prim[q] = x;
        }    
        //calculate magnitude of fgrav for tracer
        int n;
        double fg2 = 0;
        for(n=0; n<3; n++)
            fg2 += pow( c->f_grav[n], 2 );
        fgrav = sqrt( fg2 );
    }
    else{
        for( q=0; q<NUM_Q; q++ )
            tr->prim[q] = -1.0;
    }
    //Need to calculate grav_frac b4 reset v_prev
    double dvdt = get_dvdt( tr );
    if( fgrav >= 0.0 ){
        gravfrac = (dvdt - fgrav)/fgrav;
    }

    tr->gravfrac  = fabs( gravfrac );
    tr->v_prev[0] = tr->Vr;
    tr->v_prev[1] = tr->Omega*tr->R;
    tr->v_prev[2] = tr->Vz;


 */
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
    double gamma = theDomain->theParList.Adiabatic_Index; 

    struct tracer *tr = theDomain->theTracers->head;

    while( tr!=NULL ){
   	    get_tracer_cell( theDomain , tr );
	    get_local_vel( tr );
        //set_tracer_vel( tr );
        updatePhysQuants( tr, gamma );
        //calc_dPdp( theDomain );
        moveTracers( theDomain , tr , dt );
        tr = tr->next;
   }

}
