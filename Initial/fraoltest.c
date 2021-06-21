
#include "../paul.h"
static const double Gm_s = 1;

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   
   double r   = x[0];
   double phi = x[1];


   double X = 0.0;
   if( cos(phi) > 0.0 ) X = 1.0;

   double omega, disk, angular_v, radial_v;
   double w, R0, l0,rho0, rhoatm, patm, l, e, phip;
   w = 0.1;
   l0 = 0.8;
   patm = 1.0e-3;
   rhoatm = 1.0e-3;
   rho0 = 1 - rhoatm;

   e = 0.3;
   phip = 2*(22.0/7.0);
   l = r*(1+e*cos(phi-phip));
   
   R0 = l/(1+(e*cos(phi-phip)));
   angular_v = (sqrt(Gm_s*l)/(r*r));
   
   radial_v = ((r*r)/l)*(e*sin(phi-phip))*angular_v;
   
   disk = rhoatm + (rho0*exp(-((l-l0)*(l-l0))/(2*w*w)));


   omega = sqrt(Gm_s/(r*r*r));

//   omega = Om*exp(-.5*r*r/R/R);
//   Pp = P0 - .5*rho*Om*Om*R*R*exp(-r*r/R/R);

   prim[RHO] = disk;
   prim[PPP] = patm;
   prim[URR] = radial_v;
   prim[UPP] = angular_v;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
