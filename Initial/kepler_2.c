
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double q_planet = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   q_planet	= theDomain->theParList.Mass_Ratio;
}

void initial( double * prim , double * x ){

   double r = x[0];

   double rho = 1.0;
   double Pp = rho/Mach/Mach;
   double mu = q_planet/(1.+q_planet);
   double omega02 = mu/(r*r*r);
   double omegaP2 = 0.0;

   double omega = sqrt( omega02 - omegaP2 );

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0; //-1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
