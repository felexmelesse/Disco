
#include "../paul.h"

static int visc_flag = 0;
static double nu = 0.0;
static double r0 = 0.0;
static double dr = 0.0;
static int prof = 0;
static double sig0 = 0.0;
static double sig1 = 0.0;
static double mach = 0.0;
static double gam = 0.0;
static double M = 0.0;
static double a = 0.0;
static double q = 0.0;

void setICparams( struct domain * theDomain ){
    visc_flag = theDomain->theParList.visc_flag;
    nu = theDomain->theParList.viscosity;
    prof = theDomain->theParList.initPar0;
    r0 = theDomain->theParList.initPar1;
    dr = theDomain->theParList.initPar2;
    sig0 = theDomain->theParList.initPar3;
    sig1 = theDomain->theParList.initPar4;
    mach = theDomain->theParList.Disk_Mach;
    gam = theDomain->theParList.Adiabatic_Index;
    
    M = 1.0;
    q = theDomain->theParList.Mass_Ratio;
    if(strcmp(PLANETS, "cm") == 0)
        a = 1.0;
    else
        a = 0.0;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double rho, drhodr;
   double rm = r0-dr;
   double rp = r0+dr;

   
   if(prof == 1) //Gaussian
   {
       rho = sig0 + (sig1-sig0)*exp(-0.5*(r-r0)*(r-r0)/(dr*dr));
       drhodr = -(sig1-sig0)*(r-r0)/(dr*dr);
   }
   else if (prof == 2)  //Smoothed top hat, smoothing length = 0.2*r0
   {
       double d = 0.2*r0;
       rho = sig0 + (sig1-sig0) * 0.5*(1.0+tanh(-(r-rm)*(r-rp)/ (d*d)));
       double coshr = cosh((r-rm)*(r-rp)/(d*d));
       drhodr = -(sig1-sig0) * (r-0.5*(rm+rp)) / (coshr*coshr*d*d);
   }
   else
   {
       rho = sig1;
       drhodr = 0.0;
   }

   double cs20 = 1.0/r0 / (mach*mach);
   double Pp0 = sig1 * cs20 / gam;

   double Pp = Pp0 * pow(rho / sig1, gam);
   double dPdr = (gam*Pp/rho) * drhodr;

   double qp1 = 1+q;
   double y = a/((1+q)*r);

   double gr = -M/(r*r) * (1 + 3./4. * q * y*y
                        + 45./64. * q*(q*q-q+1) * y*y*y*y);
   double omega = sqrt((dPdr/rho - gr) / r);
   double vr = visc_flag ? -1.5*nu/rho/r : 0.0;

   double X = 0.0; 
   if( r>rm && r<rp ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
