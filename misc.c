#include "paul.h"
#include <string.h>


double get_dA( double * , double * , int );
double get_dV( double * , double * );
double get_moment_arm( double * , double * );

void clean_pi( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double phi_max = theDomain->phi_max;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         double phi = c->piph;
         while( phi > phi_max ){ phi -= phi_max; }
         while( phi < 0.0     ){ phi += phi_max; }
         c->piph = phi;
      }
   }
   int p;
   for( p=0 ; p<Npl ; ++p ){
      struct planet * pl = theDomain->thePlanets + p;
      double phi = pl->phi;
      while( phi > phi_max ){ phi -= phi_max; pl->RK_phi -= phi_max; }
      while( phi < 0.0     ){ phi += phi_max; pl->RK_phi += phi_max; }
      pl->phi = phi;
   }

   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      double phi = tr->Phi;
      while( phi > phi_max ){ phi -= phi_max; tr->RK_phi -= phi_max; }
      while( phi < 0.0     ){ phi += phi_max; tr->RK_phi += phi_max; }
      tr->Phi = phi;
      tr = tr->next;
   }

}

double mindt( double * , double , double * , double * );

double getmindt( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double dt = theDomain->theParList.maxDT / theDomain->theParList.CFL;
   if(dt <= 0.0)
       dt = 1.0e100; //HUGE_VAL

   int i,j,k;
   for( j=1 ; j<Nr-1 ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            if(!(c->real))
                continue;
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};
            double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
            int im = i-1;
            if( i==0 ) im = Np[jk]-1;
            double wm = theCells[jk][im].wiph;
            double wp = c->wiph;
            double w = .5*(wm+wp);
            double dt_temp = mindt( c->prim , w , xp , xm );
            //printf("jk, i, dt, dt_temp: %d, %d, %g, %g\n", jk,i,dt,dt_temp);
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   dt *= theDomain->theParList.CFL;
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );
   
   return( dt );
}

//void initial( double * , double * );
void prim2cons( double * , double * , double * , double );
void cons2prim( double * , double * , double * , double, struct planet * );
void restart( struct domain * );
/*
void clear_w( struct domain * theDomain ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         theDomain->theCells[jk][i].wiph = 0.0;
      }
   }
}*/

double get_omega( double * , double * );
double mesh_om( double );

void set_wcell( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double zm = z_kph[k-1];
         double zp = z_kph[k];
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * cL = &(theCells[jk][i ]);
            double w = 0.0;
            if( mesh_motion ){
               int ip = (i+1)%Np[jk];
               double phip = cL->piph;
               double phim = phip-cL->dphi;
               double xp[3] = {rp,phip,zp};
               double xm[3] = {rm,phim,zm};
               double r = get_moment_arm( xp , xm );
               double x[3] = {r, 0.5*(phim+phip), 0.5*(zm+zp)};
               double wL = get_omega( cL->prim , x );

               struct cell * cR = &(theCells[jk][ip]);
               phip = cR->piph;
               phim = phip-cR->dphi;
               xp[1] = phip;
               xm[1] = phim;
               r = get_moment_arm( xp , xm );
               x[0] = r;
               x[1] = 0.5*(phim+phip);
               double wR = get_omega( cR->prim , x );

               w = .5*(wL + wR);
            }
            cL->wiph = w;
         }
      }
   }
   if( mesh_motion == 3 ){
      for( j=0 ; j<Nr ; ++j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            double w = 0.0;
            for( i=0 ; i<Np[jk] ; ++i ){
               w += theCells[jk][i].wiph;
            }
            w /= (double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i ){
               theCells[jk][i].wiph = w;
            }
         }
      }
   }
   if( mesh_motion == 4 ){
      for( j=0 ; j<Nr ; ++j ){
         double r = .5*(r_jph[j]+r_jph[j-1]);
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               theCells[jk][i].wiph = r*mesh_om(r);
            }
         }
      }
   }

}

//void initial( double * , double * );
void clear_cell( struct cell * );

void regrid( struct domain * theDomain ){
/*
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double jthresh = 0.1;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
   for( k=0 ; k<Np ; ++k ){
      int jk = j+Nt*k;
      for( i=2 ; i<Nr[jk]-1 ; ++i ){
         double jump = log( theCells[jk][i-1].prim[RHO]/theCells[jk][i+1].prim[RHO] );
         if( fabs(jump) > jthresh ){
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = (c-1)->riph;
            int Nsplit = (int)(fabs(jump/jthresh));
            if(Nsplit>5) Nsplit=5;
            if(Nsplit>3) printf("r=%.2e jump=%.2e split = %d (No Cause For Alarm)\n",theCells[jk][i].riph,jump,Nsplit);
            int blocksize = (Nr[jk]-1) - i;
            Nr[jk] += Nsplit;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            memmove( theCells[jk]+i+1+Nsplit , theCells[jk]+i+1 , blocksize*sizeof(struct cell) );
            int l;
            for( l=0 ; l<Nsplit+1 ; ++l ){
               int m = l+i;
               struct cell * cnew = theCells[jk]+m;
               clear_cell( cnew );
               double rlp = rm + (rp-rm)*( (double)(l+1)/(double)(Nsplit+1) );
               double rlm = rm + (rp-rm)*( (double)(l  )/(double)(Nsplit+1) );
               double xp[3] = {rlp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rlm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( cnew->prim , x );
               cnew->riph = rlp;
               cnew->dr = rlp-rlm;
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rlp*rlp*rlp-rlm*rlm*rlm)/(rlp*rlp-rlm*rlm);
               prim2cons( cnew->prim , cnew->cons , rr , dV );
            }
            i += Nsplit;
         }
      }
   }
   }
*/
}


void adjust_RK_cons( struct domain * theDomain , double RK ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk,q;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         for( q=0 ; q<NUM_Q ; ++q ){
            c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
         }
         for( q=0 ; q<NUM_FACES ; ++q ){
            c->Phi[q] = (1.-RK)*c->Phi[q] + RK*c->RK_Phi[q];
         }
      }
   }
}

void planet_RK_adjust( struct planet * , double );

void adjust_RK_planets( struct domain * theDomain , double RK ){
   int Npl = theDomain->Npl;
   int p;
   for( p=0 ; p<Npl ; ++p ){
      planet_RK_adjust( theDomain->thePlanets+p , RK );
   }
}

void tracer_RK_adjust( struct tracer * , double);

void adjust_RK_tracers( struct domain *theDomain , double RK){

   struct tracer *tr = theDomain->theTracers->head;
   while( tr!= NULL ){
     tracer_RK_adjust( tr , RK );
     tr = tr->next;
  }
}

void move_cells( struct domain * theDomain , double dt){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         c->piph += c->wiph*dt;
      }
   }
}

double get_dp( double , double );

void calc_dp( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         int im = i-1;
         if( i == 0 ) im = Np[jk]-1;
         double phim = theCells[jk][im].piph;
         double phip = theCells[jk][i ].piph;
         double dphi = get_dp(phip,phim);
         theCells[jk][i].dphi = dphi;
      }
   }
}

void calc_prim( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      double rm = r_jph[j-1];
      double rp = r_jph[j];
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {rp,phip,z_kph[k]  };
            double xm[3] = {rm,phim,z_kph[k-1]};
            double r = get_moment_arm( xp , xm );
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), 0.5*(z_kph[k]+z_kph[k-1])};
            cons2prim( c->cons , c->prim , x , dV, theDomain->thePlanets );
         }
      }
   }
}

void calc_cons( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k;
   for( j=0 ; j<Nr ; ++j ){
      double rm = r_jph[j-1];
      double rp = r_jph[j];
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {rp,phip,z_kph[k]  };
            double xm[3] = {rm,phim,z_kph[k-1]};
            double r = get_moment_arm( xp , xm );
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), 0.5*(z_kph[k]+z_kph[k-1])};
            prim2cons( c->prim , c->cons , x , dV );
         }
      }
   }
}

void plm_phi( struct domain * );
void riemann_phi( struct cell * , struct cell * , double * , double, struct planet * );

void phi_flux( struct domain * theDomain , double dt ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int i,j,k;
   plm_phi( theDomain );
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            int ip = (i+1)%Np[jk];
            struct cell * cp = theCells[jk];
            double phi = cp[i].piph;
            double xp[3] = {r_jph[j]  ,phi,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phi,z_kph[k-1]};
            double r = get_moment_arm(xp,xm);
            double dA = get_dA(xp,xm,0);
            double x[3] = {r, phi, 0.5*(z_kph[k-1]+z_kph[k])};
            //Set cell viscosity here?
            riemann_phi( &(cp[i]) , &(cp[ip]) , x , dA*dt, theDomain->thePlanets );
         }
      }
   }

}

void buildfaces( struct domain * , int , int );
void plm_trans( struct domain * , struct face * , int , int );
void riemann_trans( struct face * , double , int, struct planet * );

void trans_flux( struct domain * theDomain , double dt , int dim ){

   int Nf;
   struct face * theFaces;
   if( dim==1 ){
      Nf = theDomain->fIndex_r[theDomain->N_ftracks_r];
      theFaces = theDomain->theFaces_1;
   }else{
      Nf = theDomain->fIndex_z[theDomain->N_ftracks_z];
      theFaces = theDomain->theFaces_2;
   }

   plm_trans( theDomain , theFaces , Nf , dim );

   int f;
   for( f=0 ; f<Nf ; ++f ){
      riemann_trans( theFaces + f , dt , dim, theDomain->thePlanets );
   }

}


void setup_faces( struct domain * theDomain , int dim ){

   struct face ** theFaces;
   int NN;
   int * nn;
   if( dim==1 ){
      theFaces = &(theDomain->theFaces_1);
      nn = theDomain->fIndex_r;
      NN = theDomain->N_ftracks_r;
   }else{
      theFaces = &(theDomain->theFaces_2);
      nn = theDomain->fIndex_z;
      NN = theDomain->N_ftracks_z;
   }

   buildfaces( theDomain , dim , 0 );
   int Nf = nn[NN];
   *theFaces = (struct face *) malloc( Nf*sizeof(struct face) );
   buildfaces( theDomain , dim , 1 );

}

void reset_fgrav( double *grav ){

    grav[0] = 0.0;
    grav[1] = 0.0;
    grav[2] = 0.0;
}

void source( struct domain *, double * , double * , double * , double * , double, double, int );
void planet_src( struct planet * , double * , double *, double * , double * , double * , double );
void omega_src( double * , double * , double * , double * , double );
void density_sink( struct domain *, double *, double *, double *, double *, double ,double );
void density_sink_yike( struct domain *, double *, double *, double *, double *, double ,double );
void damping( double *, double *, double *, double, double );

void add_source( struct domain * theDomain , double dt, int last_step ){

   struct cell ** theCells = theDomain->theCells;
   struct planet * thePlanets = theDomain->thePlanets;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;


   //Reset torques before calculating in sink fxn
   //double *torques = theDomain->torques;
   //torques[0] = 0.0;
   //torques[1] = 0.0;
   memset( theDomain->torques, 0.0, 2*sizeof(double) );

   int i,j,k,p;
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = c->piph - c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV(xp,xm);
            reset_fgrav( c->f_grav );
            source( theDomain, c->prim , c->cons , xp , xm , dV, dt, last_step );
            for( p=0 ; p<Npl ; ++p ){
               planet_src( thePlanets+p , c->prim , c->cons, c->f_grav , xp , xm , dV*dt );
            }
            omega_src( c->prim , c->cons , xp , xm , dV*dt );
            //density_sink( theDomain, c->prim, c->cons, xp, xm, dV, dt );
            density_sink_yike( theDomain, c->prim, c->cons, xp, xm, dV, dt ); 
            damping( c->cons, xp, xm, dV, dt );
         }
      }
   }

}

/*
void density_sink( struct domain *theDomain, double *prim, double *cons, double *xp, double *xm, double dV ){

   double rp = x[0];
   double rm = x[1];
   double r  = 0.5*(rp+rm);
   double dphi = get_dp(xp[1],xm[1]);
   double phi  = xm[1] + 0.5*dphi;
   double nearest_pl = nearest_planet_dist( theDomain, r, phi );
   if( nearest_pl < r_sink ){
      rho = 
   }

}
*/

void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , struct cell * sweep , int j , int k ){

   int Nr = theDomain->Nr;
   int jk = j + Nr*k;
   int Np = theDomain->Np[jk];
   double * r_jph = theDomain->r_jph;
   double dr = r_jph[j]-r_jph[j-1];
   double r  = r_jph[j];
   double Long = 0.0;
   double Short = 0.0;
   int iLong = 0;
   int iShort = 0;
   int i;
   for( i=0 ; i<Np ; ++i ){
      double dx = dr;
      double dy = r*sweep[i].dphi;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; }
      if( Short < s ){ Short = s; iShort = i; }
   }
   *L  = Long;
   *iL = iLong;
   *S  = Short;
   *iS = iShort;

}

void AMRsweep( struct domain * theDomain , struct cell ** swptr , int jk ){

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j = jk%Nr;
   int k = jk/Nr;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * sweep = *swptr;
   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , sweep , j , k );
   //printf("Long = %e, Short = %e\n",L,S);

   if( S>MaxShort ){
      int iSp = (iS+1)%Np[jk];

      //Possibly shift iS backwards by 1
      int iSm = iS-1;
      if( iSm == -1 ) iSm = Np[jk]-1;
      double dpL = sweep[iSm].dphi;
      double dpR = sweep[iSp].dphi;
      if( dpL < dpR ){
         --iS;
         --iSm;
         --iSp;
         if( iS  == -1 ) iS  = Np[jk]-1;
         if( iSm == -1 ) iSm = Np[jk]-1;
         if( iSp == -1 ) iSp = Np[jk]-1;
      }

      //Remove Zone at iS+1
      sweep[iS].dphi  += sweep[iSp].dphi;
      sweep[iS].piph   = sweep[iSp].piph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iS].cons[q]   += sweep[iSp].cons[q];
         sweep[iS].RKcons[q] += sweep[iSp].RKcons[q];
      }
      double phip = sweep[iS].piph;
      double phim = phip - sweep[iS].dphi;
      double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double r  = get_moment_arm( xp , xm );
      double dV = get_dV( xp , xm );
      double x[3] = {r, 0.5*(phim+phip), 0.5*(z_kph[k-1]+z_kph[k])};
      cons2prim( sweep[iS].cons , sweep[iS].prim , x , dV, theDomain->thePlanets );
      //Shift Memory
      int blocksize = Np[jk]-iSp-1;
      if( iSp != Np[jk]-1 ) memmove( sweep+iSp , sweep+iSp+1 , blocksize*sizeof(struct cell) );
      Np[jk] -= 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong ){
      Np[jk] += 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      int blocksize = Np[jk]-iL-1;
      memmove( sweep+iL+1 , sweep+iL , blocksize*sizeof(struct cell) );

      double dphi = sweep[iL].dphi;
      double phip = sweep[iL].piph;
      double phim = phip - dphi;
      double phi0 = .5*(phip+phim);

      sweep[iL].piph   = phi0;
      sweep[iL].dphi   = .5*dphi;
      sweep[iL+1].dphi = .5*dphi;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iL].cons[q]     *= .5;
         sweep[iL].RKcons[q]   *= .5;
         sweep[iL+1].cons[q]   *= .5;
         sweep[iL+1].RKcons[q] *= .5;
      }

      double xp[3] = {r_jph[j]  ,phi0,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double dV = get_dV( xp , xm );
      double r  = get_moment_arm( xp , xm );
      double x[3] = {r, 0.5*(xp[1]+xm[1]), 0.5*(xp[2]+xm[2])};
      cons2prim( sweep[iL].cons , sweep[iL].prim , x , dV, theDomain->thePlanets );

      xp[1] = phip;
      xm[1] = phi0;
      dV = get_dV( xp , xm );
      r  = get_moment_arm( xp , xm );
      x[0] = r;
      x[1] = 0.5*(xp[1]+xm[1]);
      cons2prim( sweep[iL+1].cons , sweep[iL+1].prim , x , dV, theDomain->thePlanets );

   }
}

void AMR( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      AMRsweep( theDomain , theCells+jk , jk );
   }
}
   
int getN0( int, int, int );

void setProcessCoords( struct domain *theDomain ){
   //TODO:Make this just a function call for tracer initialization and 
   //     parallel exchange--it's redundant save bounds

  //Sets coordinate bounds of each process without ghost cells
  //Built for use in implementing tracers
  int Nr = theDomain->Nr;
  int Nz = theDomain->Nz;
  int Ng = theDomain->Ng;
  int *dim_rank = theDomain->dim_rank;
  int *dim_size = theDomain->dim_size;
  double *r_jph = theDomain->r_jph;
  double *z_kph = theDomain->z_kph;
  int Z_Periodic = theDomain->theParList.Z_Periodic;

  //Might need to fix this to adapt to different BCs?
  double r0, rf;
  if( dim_rank[0]!=0) 
     r0 = r_jph[Ng-1];
  else
     r0 = r_jph[-1];

  if( dim_rank[0]!=dim_size[0]-1)
     rf = r_jph[(Nr-1)-Ng];
  else
     rf = r_jph[Nr-1];

  double z0, zf;
  if( dim_rank[1]!=0 || Z_Periodic)
     z0 = z_kph[Ng-1];
  else
     z0 = z_kph[-1];

  if( dim_rank[1]!=dim_size[1]-1 || Z_Periodic )
     zf = z_kph[(Nz-1)-Ng];
  else
     zf = z_kph[Nz-1];

  theDomain->r0 = r0;
  theDomain->z0 = z0;
  theDomain->delr = rf-r0;
  theDomain->delz = zf-z0;
}

void print_welcome(){
    printf("\nDisco!\n\n");
    printf("Git Version: %s\n\n", GIT_VERSION);
    printf("*Compile-time Options*\n");
    printf("HYDRO: %s\n", HYDRO);
    printf("INITIAL: %s\n", INITIAL);
    printf("BOUNDARY: %s\n", BOUNDARY);
    printf("OUTPUT: %s\n", OUTPUT);
    printf("RESTART: %s\n", RESTART);
    printf("PLANETS: %s\n", PLANETS);
    printf("HLLD: %s\n", HLLD);
    printf("ANALYSIS: %s\n", ANALYSIS);
    printf("METRIC: %s\n", METRIC);
    printf("FRAME: %s\n", FRAME);
    printf("NUM_C: %d\n", NUM_C);
    printf("NUM_N: %d\n", NUM_N);
    printf("CT_MODE: %d\n", CT_MODE);
    printf("\n");
}
 
