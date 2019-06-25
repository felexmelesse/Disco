
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double rho0 = 0.4;
   double mu = 0.1;

   double rho = (1.0/(mu*sqrt(2.0*M_PI)))*exp(-0.5*((r-rho0)*(r-rho0))/(mu*mu));
   double Pp = rho/Mach/Mach;
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = 0.0;

   double omega = sqrt( omega02 - omegaP2 );

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
