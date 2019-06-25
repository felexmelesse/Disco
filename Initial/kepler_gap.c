
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double q_planet	= 0.0;
static double a		= 1.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   q_planet	= theDomain->theParList.Mass_Ratio;
}

void initial( double * prim , double * x ){

   double r		= x[0];
   double phi		= x[1];
   double rx		= r*cos(phi);
   double ry		= r*sin(phi);
   double R		= sqrt( (a-rx)*(a-rx) + ry*ry );

   double mu		= q_planet/(1.+q_planet);
   double omega02	= mu/(r*r*r);
   double omegaP2	= 0.0;

   double omega		= sqrt( omega02 - omegaP2 );

   double alpha_visc	= nu*Mach*Mach/(a*a*omega);
   double K		= q_planet*q_planet*Mach*Mach*Mach*Mach*Mach/alpha_visc;
   
   double f_0		= 0.45;
   double tau_sh	= 1.89 + 0.53/(q_planet*Mach*Mach*Mach);
   double tau_r		= 0.3363585661*pow(fabs(1.5*Mach*((R/a) - 1.0)), 2.5);
   double f_r;
   if (tau_r < tau_sh){
	f_r	= f_0;
   }else{
	f_r	= f_0*sqrt(tau_sh / tau_r);
   }
   
   double rho_0		= 1.0;
   double nom		= f_r*K/(3.0*M_PI);
   double denom		= 1.0 + f_0*K/(3.0*M_PI);
   double rho		= rho_0 * (1.0 - (nom/denom)*sqrt(a/R));
   double Pp		= rho/Mach/Mach;

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0; //-1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
