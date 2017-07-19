
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(8);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* Generic Diagnostics for Euler*/

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];
   double vz = prim[UZZ];
   double Pp = prim[UPP];

   Qrz[0] = rho;
   Qrz[1] = 2.*M_PI*r*rho*vr;
   Qrz[2] = 2.*M_PI*r*rho*vz;
   double Fr,Fp,Fz;
   Fp = 0.0;
   if( theDomain->Npl > 1 ){
      struct planet * pl = theDomain->thePlanets+1;
      planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 0 );
   }
   Qrz[3] = 2.*M_PI*r*rho*(r*Fp);
   Qrz[4] = 2.*M_PI*r*rho*( r*r*omega*vr );
   Qrz[5] = omega;
   Qrz[6] = 2*M_PI*r * rho*r*r*omega;
   Qrz[7] = Pp;
}

