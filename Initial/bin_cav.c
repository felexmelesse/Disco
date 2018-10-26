
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double rho_fl = 0.0;
static double viscosity = 0.0;
static int alpha = 0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   rho_fl = theDomain->theParList.Density_Floor;
   viscosity = theDomain->theParList.viscosity;
   alpha = theDomain->theParList.alpha_flag;
}


double get_cs2( double, double );

void initial( double * prim , double * x ){
   
   double r = x[0];
   double z = x[2];
   double phi = x[1];

   double d    = 10.0;
   double r0   = 2.5; 
   double sint = z/sqrt(r*r+z*z);
   double cs2  = get_cs2( r, phi );

   double rho, Pp;
   //rho = 1.0*exp(-sint*sint*Mach*Mach);
   rho = pow(r,-0.5)*exp( -pow(r/r0,-d) ) + rho_fl;
   
   Pp  = rho*cs2/gam;
   //Pp  = 0.01*pow(r,-3./2.)*exp( -pow(r/r0,-d) );
   
   //double n = 1.5;
   //double omega = ( pow( r , n-1.5 ) + 1. )/( pow( r , n ) + 1. );
   //double omega = 1./pow( pow( r , 1.5*n ) + 0.3 , 1./n );
   double omega2  = pow( r , -3. );
   double omega2P = 0.0*cs2/r/r;
   double omega = sqrt( omega2 - omega2P );

   if( omega > 1. ) omega = 1; 

   double X = 0.0; 
   if( r > 1.0 ) X = 1.0; 

   double nu = viscosity;
   if( alpha ){
     double alpha = viscosity;
     nu = alpha*cs2*pow(r,1.5);
   }

   double vr0 = -1.5*nu/r; 
   double vr  = vr0*exp( -pow(r/r0,-d) );

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;   
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
  
}
