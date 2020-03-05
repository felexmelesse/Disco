
#include "../paul.h"

static double nu = 0.0;

void setICparams( struct domain * theDomain ){
   nu = theDomain->theParList.viscosity;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double c = 0.0;

   double rho = 1. + 1./sqrt(r) + c*exp(-200.*pow(r-.5,2.));
   double Pp  = 0.0001;
   double omega = 1.0/pow(r,1.5);

   double drhodr = -0.5*pow(r, -1.5) - c*400*(r-0.5)*exp(-200*pow(r-.5,2));

   double vr = -1.5*nu/r * (1 + 2*r*drhodr/rho);

   double X = 0.0; 
   if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
