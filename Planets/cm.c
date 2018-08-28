
#include "../paul.h"

static double q_planet = 1.0;
static double smooth = 0.0;
static double Mach = 1.0;
static double tt = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   smooth = 0.5/Mach;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
 
   double q = q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = 1.-mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = om; 
   thePlanets[0].r     = a*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].eps   = smooth;// + 1.0;//0.025; 

   thePlanets[1].M     = mu;  
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om;  
   thePlanets[1].r     = a*(1.-mu); 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = smooth;// + 1.0;//0.025;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;
   double eps = smooth + exp(-t/30.);
   thePlanets[0].eps = eps;
   thePlanets[1].eps = eps;

   tt = t;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

double phigrav( double, double, double );
double get_planet_2dist( double, double, double, double );

double get_bin_potential( double r, double phi ){
   double a  = 1.0; 
   double q = q_planet;
   double om = pow( a, -1.5 );
   double mu = q/(1.+q);
   
   double M1 = 1.-mu;
   double r1 = a*mu;
   double phi1 = om*tt - M_PI;
   double d1 = get_planet_2dist( r, phi, r1, phi1 );

   double M2 = mu;
   double r2 = a*(1.-mu);
   double phi2 = om*tt;
   double d2 = get_planet_2dist( r, phi, r2, phi2 );

   double cut = 0.05;
   if( d2 < cut ){
       d2 = cut;
       //M2 = 1.0; //This is just to match Yike's....
   }
   if( d1 < cut ){
       d1 = cut;
       //M1 = 1.0;
   }

   double pot1 = M1/d1;    //phigrav( M1, d1, 0.0 ); //no smoothing
   double pot2 = M2/d2;    //phigrav( M2, d2, 0.0 );
   //printf("Time: %f\n", tt );
   //printf("pot1: %f, pot2: %f\n", pot1, pot2 );
   return pot1 + pot2;
}

