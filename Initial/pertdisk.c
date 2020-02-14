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

   double rho = 1.0 + 0.2*exp(-(r-10.0)*(r-10.0));
   double Pp  = rho/(Mach*r);
   double omega = pow(r,-1.5);

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;

}
