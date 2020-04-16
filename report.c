#include "paul.h"
#include "geometry.h"

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

void report( struct domain * theDomain ){

   double t = theDomain->t;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int rank = theDomain->rank;
#if USE_MPI
   MPI_Comm grid_comm = theDomain->theComm;
#endif

   double gamma_law = theDomain->theParList.Adiabatic_Index;

   struct planet * thePlanets = theDomain->thePlanets;
   int Npl = theDomain->Npl;

   //double r_p = 0.0;
   //if( Npl > 1 ) r_p = thePlanets[1].r;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int jmin = NgRa;
   int jmax = Nr-NgRb;
   int kmin = NgZa;
   int kmax = Nz-NgZb;

   int j,k,i;
   double L1_isen = 0.0;
   double L1_rho  = 0.0;
   double L1_P    = 0.0;
   double L1_B    = 0.0;
   double Br2     = 0.0;
   double B2      = 0.0;
   double Power  = 0.0;
   double Torque = 0.0;
   double Torque2 = 0.0;
   double Fr=0.0;
   double PsiR = 0.0;
   double PsiI = 0.0;
   double Vol = 0.0;
   double rho_min = HUGE_VAL;
   double rhoavg_min = HUGE_VAL;
   double Mass = 0.0;
   double Mdot = 0.0;

   double BrBp = 0.0;
   double PdV  = 0.0;

   double * M_acc, * L_pls, * Ls_pls, *kin_pls, *therm_pls, *Lt_pls;
   M_acc = calloc(Npl, sizeof(double) );
   L_pls = calloc(Npl, sizeof(double) );
   Ls_pls = calloc(Npl, sizeof(double) );
   Lt_pls = calloc(Npl, sizeof(double) );
   kin_pls = calloc(Npl, sizeof(double) );
   therm_pls = calloc(Npl, sizeof(double) );

   double S_R = 0.0;
   double S_0 = 0.0;

   //double T_cut[10];
   //double P_cut[10];
   //for( j=0 ; j<10 ; ++j ){ T_cut[j]=0.;  P_cut[j]=0.; }

   for( j=0; j<Npl; ++j){
      //M_acc[j] = 0.5*thePlanets[j].dM + 0.5*thePlanets[j].RK_dM;
      //L_pls[j] = 0.5*thePlanets[j].L + 0.5*thePlanets[j].RK_L;
      //Ls_pls[j] = 0.5*thePlanets[j].Ls + 0.5*thePlanets[j].RK_Ls;
      //therm_pls[j] = 0.5*thePlanets[j].therm + 0.5*thePlanets[j].RK_therm;
      //kin_pls[j] = 0.5*thePlanets[j].kin + 0.5*thePlanets[j].RK_kin;

      M_acc[j] = thePlanets[j].dM;
      L_pls[j] = thePlanets[j].L;
      Ls_pls[j] = thePlanets[j].Ls;
      therm_pls[j] = thePlanets[j].therm;
      kin_pls[j] = thePlanets[j].kin;
      Lt_pls[j] = thePlanets[j].Ltorque;

      thePlanets[j].dM = 0.0;
      thePlanets[j].RK_dM = 0.0;
      thePlanets[j].L = 0.0;
      thePlanets[j].RK_L = 0.0;
      thePlanets[j].Ls = 0.0;
      thePlanets[j].RK_Ls = 0.0;
      thePlanets[j].therm = 0.0;
      thePlanets[j].RK_therm = 0.0;
      thePlanets[j].kin = 0.0;
      thePlanets[j].RK_kin = 0.0;
      thePlanets[j].Ltorque = 0.0;
      thePlanets[j].RK_Ltorque = 0.0;
  }
   for( j=jmin ; j<jmax ; ++j ){
      double rho0 = 1.0;//pow( r , -1.5 );
      double rho_avg = 0.0;
      double Vol_avg = 0.0;
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phi = c->piph - .5*c->dphi;
            double Pp  = c->prim[PPP];
            double rho = c->prim[RHO];

            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double r = get_centroid( xp[0] , xm[0], 1 );

            PsiR += rho*dV*cos(phi);
            PsiI += rho*dV*sin(phi);
 
            L1_isen += fabs(Pp/pow(rho,gamma_law)-1.)*dV;
            L1_rho  += fabs(rho/rho0-1.)*dV;
            L1_P    += fabs(Pp/pow(rho,5./3.)/0.01-1.)*dV;
            if( NUM_C > BZZ ){
               double Br = c->prim[BRR];
               double Bp = c->prim[BPP];
               double Bz = c->prim[BZZ];
               L1_B += fabs(Br)*dV;
               Br2 += .5*Br*Br*dV;
               BrBp += Br*Bp*dV;
               B2  += .5*(Br*Br+Bp*Bp+Bz*Bz)*dV;
            }
            PdV += Pp*dV;
            Mdot += 2.*M_PI*r*rho*dV*c->prim[URR];
            Vol += dV;

            if( rho_min > rho/rho0 ) rho_min = rho/rho0;
            rho_avg += rho*dV;
            Vol_avg += dV;
            Mass += rho*dV;

            S_R += pow(rho,4.)*r*dV;
            S_0 += pow(rho,4.)*dV;

            if( Npl > 1 ){
               double fr,fp,fz,fp2;
               double rp = thePlanets[1].r;
               double om = thePlanets[1].omega;
               double vr = thePlanets[1].vr;
               //double mp = thePlanets[1].M;
               planetaryForce( thePlanets   , r , phi , 0.0 , &fr , &fp2 , &fz , 1 );
               planetaryForce( thePlanets+1 , r , phi , 0.0 , &fr , &fp  , &fz , 1 );

               //double pp = thePlanets[1].phi;
               //double cosp = cos(phi);
               //double sinp = sin(phi);
               //double dx = r*cosp-rp*cos(pp);
               //double dy = r*sinp-rp*sin(pp);
               //double script_r = sqrt(dx*dx+dy*dy);
               //double rH = pow( thePlanets[1].M/3. , 1./3. );

               Torque -= (rho-1.0)*rp*fp*dV;
               Power  -= (rho-1.0)*( rp*om*fp + vr*fr )*dV;
               Torque2 -= (rho-1.0)*rp*fp2*dV;
               Fr -= (rho-1.0)*fr*dV;

/*
               int n_cut;
               for( n_cut=0 ; n_cut<10 ; ++n_cut ){
                  double r_cut = 0.3*((double)(n_cut+1.)/10.);
                  double scriptr2 = rp*rp + r*r - 2.*r*rp*cos(phi-pp);
                  if( scriptr2 > r_cut*r_cut ){
                     T_cut[n_cut] -= (rho-1.0)*rp*fp*dV;
                     P_cut[n_cut] -= (rho-1.0)*( rp*om*fp + vr*fr )*dV;
                  }
               }
*/
            }
         }
      }
      rho_avg /= Vol_avg;
      if( rhoavg_min > rho_avg/rho0 ) rhoavg_min = rho_avg/rho0;
   }

#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , &L1_isen , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &L1_rho  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &L1_P    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &L1_B    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Br2     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &B2      , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &BrBp    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PdV     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Vol     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Torque  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Torque2 , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Power   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Fr      , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PsiR    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &PsiI    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rho_min    , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &rhoavg_min , 1 , MPI_DOUBLE , MPI_MIN , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Mass    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &S_R     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &S_0     , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , &Mdot    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   //MPI_Allreduce( MPI_IN_PLACE , &M_acc   , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   for( j=0; j<Npl; ++j){
      MPI_Allreduce( MPI_IN_PLACE , &M_acc[j]  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
      MPI_Allreduce( MPI_IN_PLACE , &L_pls[j]  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
      MPI_Allreduce( MPI_IN_PLACE , &Ls_pls[j]  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
      MPI_Allreduce( MPI_IN_PLACE , &kin_pls[j]  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
      MPI_Allreduce( MPI_IN_PLACE , &therm_pls[j]  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   }

//   MPI_Allreduce( MPI_IN_PLACE , T_cut  , 10 , MPI_DOUBLE , MPI_SUM , grid_comm );
//   MPI_Allreduce( MPI_IN_PLACE , P_cut  , 10 , MPI_DOUBLE , MPI_SUM , grid_comm );
#endif

   L1_isen /= Vol;
   L1_rho  /= Vol;
   L1_P    /= Vol;
   L1_B    /= Vol;
   Mdot /= Vol;
   S_R /= S_0;

   double aM = BrBp/PdV;
   double bM = PdV/B2;

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      //fprintf(rFile,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
      //          t,Torque,Power,Fr,rho_min,rhoavg_min,PsiR,PsiI,Mass,Mdot,S_R,
      //          L1_rho,L1_isen,L1_B,Br2,aM,bM,M_acc);
      //fprintf(rFile,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le ",
      //          t,Torque,Power,Fr,rho_min,rhoavg_min,PsiR,PsiI,Mass,Mdot,S_R,
      //          L1_rho,L1_isen,L1_B,Br2,aM,bM);
      //fprintf(rFile,"%le %le %le ",  t,Torque,Torque2);
      fprintf(rFile,"%le ",  t);
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", Lt_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", M_acc[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", L_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", Ls_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", kin_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%le ", therm_pls[j]);
      }
      fprintf(rFile,"\n");

      //fprintf(rFile,"%e %e %e ",t,Torque,Power);
      //for( j=0 ; j<10 ; ++j ) fprintf(rFile,"%e %e ",T_cut[j],P_cut[j]);
      //fprintf(rFile,"\n");
      fclose(rFile);
   }
   free(M_acc);
}
