
#include "paul.h"

double get_signed_dp( double , double );

double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

void plm_phi( struct domain * theDomain ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double PLM = theDomain->theParList.PLM;
   int i,j,k,q;
   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            int im = i-1;
            if( im == -1 ) im = Np[jk]-1;
            int ip = (i+1)%Np[jk];
            struct cell * c  = &(theCells[jk][i]);
            struct cell * cL = &(theCells[jk][im]);
            struct cell * cR = &(theCells[jk][ip]);
            double dpL = cL->dphi;
            double dpC = c->dphi;
            double dpR = cR->dphi;
            for( q=0 ; q<NUM_Q ; ++q ){
               double pL = cL->prim[q];
               double pC = c->prim[q];
               double pR = cR->prim[q];
               double sL = pC - pL;
               sL /= .5*( dpC + dpL );
               double sR = pR - pC;
               sR /= .5*( dpR + dpC );
               double sC = pR - pL;
               sC /= .5*( dpL + dpR ) + dpC;
               c->gradp[q] = minmod( PLM*sL , sC , PLM*sR );
            }
         }
      }
   }
}

double get_dA( double * , double * , int );
double get_centroid(double , double , int );
void geom_grad(double *prim, double *grad, double *xp, double *xm, 
                double PLM, int dim, int LR);
void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, int kmin,
                        int kmax, int dim, int LR);

void plm_trans( struct domain * theDomain , struct face * theFaces , int Nf , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   double PLM = theDomain->theParList.PLM;
   
   int i,j,k,q;

   //Clear gradients
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            memset( theCells[jk][i].grad , 0 , NUM_Q*sizeof(double) );
         }
      }
   }

   if(PLM == 0.0)
       return;

   //Add weighted slopes
   int n;
   for( n=0 ; n<Nf ; ++n ){
      struct face * f  = &( theFaces[n] );
      double phi = f->cm[1];
      struct cell * cL = f->L;
      struct cell * cR = f->R;
      double dxL = f->dxL;
      double dxR = f->dxR;
      double phiL = cL->piph - .5*cL->dphi;
      double phiR = cR->piph - .5*cR->dphi;
      double dpL = get_signed_dp(phi,phiL);
      double dpR = get_signed_dp(phiR,phi);
      double dA  = f->dA;
      for( q=0 ; q<NUM_Q ; ++q ){
         double WL = cL->prim[q] + dpL*cL->gradp[q];
         double WR = cR->prim[q] - dpR*cR->gradp[q];

         double S = (WR-WL)/(dxR+dxL);
         cL->grad[q] += S*dA;
         cR->grad[q] += S*dA;
      }
   }

   //Divide by total weight
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            if( dim==1 ) xm[0] = r_jph[j];
            else xm[2] = z_kph[k];
            double dAp = get_dA(xp,xm,dim);
            if( dim==1 ){
               xp[0] = r_jph[j-1];
               xm[0] = r_jph[j-1];
            }else{
               xp[2] = z_kph[k-1];
               xm[2] = z_kph[k-1];
            }
            double dAm = get_dA(xp,xm,dim);
            double dAtot = dAp+dAm;
            if( (dim==1 && j==0   ) || (dim==2 && k==0   ) ) dAtot = dAp;
            if( (dim==1 && j==Nr-1) || (dim==2 && k==Nz-1) ) dAtot = dAm;
            for( q=0 ; q<NUM_Q ; ++q ){
               c->grad[q] /= dAtot;
            }
         }    
      }    
   }

   //Slope Limiting
   for( n=0 ; n<Nf ; ++n ){
      struct face * f  = &( theFaces[n] );
      double phi = f->cm[1];
      struct cell * cL = f->L;
      struct cell * cR = f->R;
      double dxL = f->dxL;
      double dxR = f->dxR;
      double phiL = cL->piph - .5*cL->dphi;
      double phiR = cR->piph - .5*cR->dphi;
      double dpL = get_signed_dp(phi,phiL);
      double dpR = get_signed_dp(phiR,phi);
      for( q=0 ; q<NUM_Q ; ++q ){
         double WL = cL->prim[q] + dpL*cL->gradp[q];
         double WR = cR->prim[q] - dpR*cR->gradp[q];

         double S = (WR-WL)/(dxR+dxL);
         double SL = cL->grad[q];
         double SR = cR->grad[q];
         if( S*SL < 0.0 ){
            cL->grad[q] = 0.0; 
         }else if( fabs(PLM*S) < fabs(SL) ){
            cL->grad[q] = PLM*S;
         }
         if( S*SR < 0.0 ){
            cR->grad[q] = 0.0; 
         }else if( fabs(PLM*S) < fabs(SR) ){
            cR->grad[q] = PLM*S;
         }    
      }    
   }

   // Geometric boundaries don't have ghost zones, so the gradients there need
   // to be fixed.  
   
   if(dim == 1 && theDomain->NgRa == 0)
      plm_geom_boundary(theDomain, 0, 0, 0, Nz-1, dim, 0);

   if(dim == 1 && theDomain->NgRb == 0)
      plm_geom_boundary(theDomain, Nr-1, Nr-1, 0, Nz-1, dim, 1);

   if(dim == 2 && theDomain->NgZa == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, 0, 0, dim, 0);

   if(dim == 2 && theDomain->NgZb == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, Nz-1, Nz-1, dim, 1);
}

void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, 
                        int kmin, int kmax, int dim, int LR)
{
    struct cell ** theCells = theDomain->theCells;
    int Nr = theDomain->Nr;
    int * Np = theDomain->Np;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;
    double PLM = theDomain->theParList.PLM;

    int i,j,k;

    double xp[3], xm[3];
    for(k=kmin; k<=kmax; k++)
    {
       xp[2] = z_kph[k];
       xm[2] = z_kph[k-1];
       for(j=jmin; j<=jmax; j++)
       {
           int jk = j+Nr*k;
           xp[0] = r_jph[j];
           xm[0] = r_jph[j-1];
           for(i=0; i < Np[jk]; i++)
           {
               struct cell * c = &(theCells[jk][i]);
               xp[1] = c->piph;
               xm[1] = c->piph - c->dphi;

               geom_grad(c->prim, c->grad, xp, xm, PLM, dim, LR); 
           }
       }
    }
}
