// Fishbone & Moncrief 1976 Astrophysical Journal 207, 962-976
// http://adsabs.harvard.edu/abs/1976ApJ...207..962F

#include "../paul.h"
#include "../Hydro/frame.h"

static double M = 0.0;
static double a = 0.0;
static double gam = 0.0;
static double Rin = 0.0;
static double Rmax = 0.0;
static double rho_atm = 0.0;
static double B0 = 0.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   //M = theDomain->theParList.metricPar2;
   M = 1.0;
   a = theDomain->theParList.metricPar3;
   Rin = theDomain->theParList.initPar1;
   Rmax = theDomain->theParList.initPar2;
   rho_atm = theDomain->theParList.initPar3;
   B0 = theDomain->theParList.initPar4;
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];
   double z  = x[2];
   double R = sqrt(r*r + z*z);

   double k;
   if(Rmax > 2*Rin)
       Rmax = 2*Rin;
   if(Rmax < Rin)
       Rmax = Rin;
   k = Rmax/Rin;

   double h = M/R - 0.5*M*Rmax/(r*r) - 0.5*M*(2-k)/Rin;
   double hmax = 0.5*M*(k-1)*(k-1)/Rmax;
   double rho = pow(h/hmax, 1.0/(gam-1.0));
   double l = sqrt(M*Rmax);
   double om = l/(r*r);
   double q = 1.0;

   if(h < 0.0 || rho < rho_atm)
   {
       h = M/R;
       rho = rho_atm * pow(Rin/R, 1.0/(gam-1.0));
       q = 0.0;
       om = 0.0;
   }

   double P = (gam-1)/gam * rho*h;


   prim[RHO] = rho;
   prim[PPP] = P;
   prim[URR] = 0.0;
   prim[UPP] = om;
   prim[UZZ] = 0.0;

   if(NUM_C > BZZ)
   {
       double Bz = B0;
        prim[BRR] = 0.0;
        prim[BPP] = 0.0;
        prim[BZZ] = Bz;
   }

   if( NUM_N>0 ) prim[NUM_C] = q;

}

