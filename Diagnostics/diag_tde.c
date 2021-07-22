#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(10);
}

//void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );
double phigrav( double M , double r , double eps , int type);

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
   double vp = r*omega;
   double Pp = prim[PPP];

   Qrz[0] = rho;
   Qrz[1] = rho*vr;
   Qrz[2] = rho*vp;
   Qrz[3] = rho*vp*vr;
   Qrz[4] = 0.5*rho*(vp*vp+vr*vr);
   Qrz[5] = Pp/(gamma_law-1.0);

   //eccentricity stuff
   double gx = r*cos(phi);
   double gy = r*sin(phi);
   double vx = vr*cos(phi) - vp*sin(phi);
   double vy = vr*sin(phi) + vp*cos(phi);
   double vm2 =  vx*vx + vy*vy;
   double rdv = vx*gx + vy*gy;
   Qrz[6] = (vm2 - 1./r)*gx - rdv*vx;	//e_x
   Qrz[6] *= rho;
   Qrz[7] = (vm2 - 1./r)*gy - rdv*vy;	//e_y
   Qrz[7] *= rho;

   struct planet * pl = theDomain->thePlanets;
   double ebind = 0.5*(vr*vr + vp*vp) - phigrav( pl.M , r , 0.0 , pl.type);
   double e2 = 1.0 + 2*ebind*vp*vp;
   if (e2 < 0) e2 = 0;
   Qrz[8] = rho*sqrt(e2);

   ebind = 0.5*(vr*vr + vp*vp) - 1/r;
   e2 = 1.0 + 2*ebind*vp*vp;
   if (e2 < 0) e2 = 0;
   Qrz[9] = rho*sqrt(e2);

}

