#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int alpha_flag = 0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   alpha_flag = theDomain->theParList.alpha_flag;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double R = sqrt(r*r + eps*eps);

   double omega02 = 1.0/pow(R,3.);
   double omegaP2 = (1.5/(Mach*Mach*R*R*R));

   double omega2 = fmax( (omega02 - omegaP2), 0.0 );
   double omega = sqrt(omega2);
   double cs2 = 1.0/(R*Mach*Mach);
   double visc = nu;
   if (alpha_flag == 1) visc = nu*cs2/omega;
   double rho = 1.0;
   if (nu > 0) rho = 1.0/visc;
   double Pp = rho*cs2;

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*visc/(r);
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
