
#include "paul.h"

static int meshOmChoice = 0;
static double meshOmPar = 0.0;
static int enOmChoice = 0;
static double enOmPar = 0.0;
static int cs2Choice = 0;
static double cs2Par = 0.0;
static double Omega0 = 0.0;

static double Mach = 0.0;
static double r0 = 0.0;
static double r1 = 0.0;
static double r2 = 0.0;
static double H0 = 0.0;
static double M = 0.0;

static struct planet *thePlanets = NULL;


void setOmegaParams( struct domain * theDomain ){
   meshOmChoice = theDomain->theParList.Exact_Mesh_Omega;
   meshOmPar    = theDomain->theParList.Exact_Mesh_Omega_Par;
   enOmChoice = theDomain->theParList.Energy_Omega;
   enOmPar = theDomain->theParList.Energy_Omega_Par;
   cs2Choice = theDomain->theParList.Cs2_Profile;
   cs2Par = theDomain->theParList.Cs2_Par;
   Omega0 = theDomain->theParList.RotOmega;

   Mach = theDomain->theParList.Disk_Mach;
   r0 = theDomain->theParList.initPar1; // Inner edge
   r1 = theDomain->theParList.initPar2; // Fiducial radius
   r2 = theDomain->theParList.initPar3; // Outer Edge
   H0 = theDomain->theParList.initPar4; // Scale Height
   M = theDomain->theParList.metricPar2;

   thePlanets = theDomain->thePlanets;

   if(strcmp(PLANETS, "bin_rot") == 0)
   {
       double q = theDomain->theParList.Mass_Ratio;
       M = q*1.0/(1.0+q);
   }
}

double mesh_om( double *x)
{
    double r = x[0];
    double omega;
    if(meshOmChoice == 1)
        omega = 1.0;
    else if(meshOmChoice == 2)
        omega = pow(r,-1.5);
    else if(meshOmChoice == 3)
    {
        double n = 8.0;
        omega = 1./pow( pow( r , 1.5*n ) + 1. , 1./n );
    }
    else
        omega = 0.0;
        
   return( omega );
}

double get_om( double *x ){
    double r = x[0];
    double om;

    if(enOmChoice == 1)
        om = 1.0-Omega0;

    else if(enOmChoice == 2)
        om = 1.0/pow(r,1.5)-Omega0;

    else if(enOmChoice == 3)
    {
        double n = 8.0;
        om = 1./pow( pow( r , 1.5*n ) + 1. , 1./n )-Omega0;
    }

    else if(enOmChoice == 4)
        om =  10.*exp(-.5*r*r)-Omega0;

    else
        om = 0.0-Omega0;

    return om;
}
  
double get_om1( double *x){
    double r = x[0];
    double om1;

    if(enOmChoice == 1)
        om1 = 0.0;

    else if(enOmChoice == 2)
        om1 = -1.5/pow(r,2.5);

    else if(enOmChoice == 3)
    {
        double n = 8.0;
        om1 = -1.5 * pow(r,-1+1.5*n) / pow( pow(r,1.5*n) + 1. , 1.0+1./n );
    }

    else if(enOmChoice == 4)
        om1 =  10.*r*exp(-.5*r*r);

    else
        om1 = 0.0;

    return om1;
}
  
double get_om2( double *x){
    return 0.0;
}

double get_cs2( double *x ){
    double r = x[0];
    double cs2;

    if(cs2Choice == 1)
        cs2 = 1./(Mach*Mach);
    else if(cs2Choice == 2)
    {
        double nu = .5;
        cs2 = .5/Mach/Mach/pow(r,2.*nu);
    }
    else if(cs2Choice == 3)
    {
        cs2 = M*H0*H0 / (2*r1*r1*r1);
    }
    else if(cs2Choice == 4)
    {
        double v2 = M/r;
        cs2 = v2/(Mach*Mach);
    }
    else if(cs2Choice == 5)
    {
      double r = 0.5*x[0];
      double phi = 0.5*x[1];

      double cosp = cos(phi);
      double sinp = sin(phi);
      double gx = r*cosp;
      double gy = r*sinp;
      int pi;
      double px, py, pr;
      double n = 2.0;
      double phip = 0.0;
 
      for (pi = 0; pi<2; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;
        pr = (px-gx)*(px-gx) + (py-gy)*(py-gy);

        phip += M/pow( pr + pow(thePlanets[pi].eps,n) , 1./n );
      }
      cs2 = phip/(Mach*Mach);        
    }
    else
        cs2 = 1.0;

    return cs2;
}

