
#include "../paul.h"
static const double Gm_s = 1;
static double r_out = 1;

void setICparams( struct domain * theDomain ){
   r_out = theDomain->theParList.rmax;
}

void initial( double * prim , double * x ){
   
   double r   = x[0];
   double phi = x[1];


   double X = 0.0;
//if( cos(phi) > 0.0 ) X = 1.0;

   double omega, angular_v, radial_v, E, J, rho, p, vr, E_mb, a_mb ,e_mb, M;
   double w, R0, rho0, rhoatm, patm, l, e, phip, phi0, rp, a, M_star, R_star, beta, R_T;
   //E = -0.02;
   //J = 3;
   M_star = 1.0e-6;
   R_star = 0.47;
   beta = 4;
   M = 1;
   //rp = 4;
   //ra = 20;
   w = 0.0025;
   //l0 = 0.8;
   patm = 1.0e-5;
   rhoatm = 1.0e-3;
   rho0 = 1 - rhoatm;

   R_T = R_star*pow(M/M_star,1.0/3.0);
   rp = R_T/beta;
   E_mb = -M*R_star/(R_T*R_T);
   a_mb = (R_T*R_T)/2*R_star;
   e_mb = 1-(rp/a_mb);
   a = a_mb;
   e = e_mb;
   l = a*(1-(e*e));
   phip = M_PI;
   

   J = sqrt(Gm_s*l);
   E = -0.5*(1-(e*e))*(Gm_s/l);
   
   R0 = l/(1+(e*cos(phi-phip)));
   phi0 = phip - acos(((l/r)-1)/e); // Multi by (-ve) to make sure the gas enters the grid instaed of exit
   if(phi0 - phi < -M_PI)
       phi0 += 2*M_PI;
   if(phi0-phi > M_PI)
      phi0 -= 2*M_PI;
   angular_v = (sqrt(Gm_s*l)/(r*r));
   
   radial_v = ((r*r)/l)*(e*sin(phi-phip))*angular_v;
   
//   disk = rhoatm + (rho0*exp(-((l-l0)*(l-l0))/(2*w*w)));
   //printf("phip %f, phi %f, phi0 %f\n",phip, phi, phi0);
   printf("%lf\n",a_mb);
   
   if ((r > r_out) && (phi0 - 5*w < phi) && (phi < phi0 + 5*w)){
// In the stream
      rho = 1;
      p = 1.0e-5;
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
