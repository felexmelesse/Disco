
#include "../paul.h"
static const double Gm_s = 1;

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   
   double r   = x[0];
   double phi = x[1];


   double X = 0.0;
//if( cos(phi) > 0.0 ) X = 1.0;

   double omega, disk, angular_v, radial_v, E, J, rho, p, vr, r_out;
   double w, R0, l0,rho0, rhoatm, patm, l, e, phip, phi0;
   E = -0.02;
   J = 3;
   r_out =10;
   w = 0.02;
   l0 = 0.8;
   patm = 1.0e-3;
   rhoatm = 1.0e-3;
   rho0 = 1 - rhoatm;

   l = (J*J)/(Gm_s);
   e = sqrt(1+((2*J*J*E)/(Gm_s*Gm_s)));

   //e = 0.5;
   phip = M_PI;
   //l = r*(1+e*cos(phi-phip));
   
   R0 = l/(1+(e*cos(phi-phip)));
   phi0 = phip - acos(((l/r)-1)/e); // Multi by (-ve) to make sure the gas enters the grid instaed of exit
   if(phi0 - phi < -M_PI)
       phi0 += 2*M_PI;
   if(phi0-phi > M_PI)
      phi0 -= 2*M_PI;
   angular_v = (sqrt(Gm_s*l)/(r*r));
   
   radial_v = ((r*r)/l)*(e*sin(phi-phip))*angular_v;
   
   disk = rhoatm + (rho0*exp(-((l-l0)*(l-l0))/(2*w*w)));
   //printf("phip %f, phi %f, phi0 %f\n",phip, phi, phi0);
   
   if ((r > r_out) && (phi0 - 5*w < phi) && (phi < phi0 + 5*w)){
// In the stream
      rho = 1;
      p = 0.01;
      vr = radial_v;
      omega = angular_v;
      X = 1;
   }  
   else{
// In the atm 
      rho = rhoatm;
      p = patm;
      vr = 0;   
      omega = sqrt(Gm_s/(r*r*r));  
      X = 0;
   }
   

//   omega = Om*exp(-.5*r*r/R/R);
//   Pp = P0 - .5*rho*Om*Om*R*R*exp(-r*r/R/R);

   prim[RHO] = rho;
   prim[PPP] = p;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
