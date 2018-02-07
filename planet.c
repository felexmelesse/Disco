#include "paul.h"

double PHI_ORDER = 2.0;
static double R_SCHWZ = 0.1/3.0;
static int grav2D = 0;
static int pw_flag = 1.0;

double get_dp( double , double );

void setGravParams( struct domain * theDomain ){

   grav2D = theDomain->theParList.grav2D; 

}

double phigrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M/pow( pow(r,n) + pow(eps,n) , 1./n ) ) ;
}

double fgrav( double M , double r , double eps ){
   double n = PHI_ORDER;
   return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n),1.+1./n) );
}

double phigrav_pw( double M, double r ){
    return( M/(r-R_SCHWZ) );
}

double fgrav_pw( double M, double r ){
    return( M/pow(r-R_SCHWZ, 2) );
}

void adjust_gas( struct planet * pl , double * x , double * prim , double gam ){

   double r   = x[0];
   double phi = x[1];

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy);

   double pot = 0.0;
   double r_cut = 1.5*R_SCHWZ;
   if( pw_flag ){
        if( script_r < r_cut )
            script_r = r_cut;
        pot = phigrav_pw( pl->M, script_r );
   }
   else{
        pot = phigrav( pl->M , script_r , pl->eps );
   }

   double c2 = gam*prim[PPP]/prim[RHO];
   double factor = 1. + (gam-1.)*pot/c2;

   prim[RHO] *= factor;
   prim[PPP] *= pow( factor , gam );

}

void planetaryForce( struct planet * pl , double r , double phi , double z , double * fr , double * fp , double * fz , int mode ){

   if(grav2D)
      z = 0.0;

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy+z*z);
   double script_r_perp = sqrt(dx*dx+dy*dy);

   double r_cut = 1.5*R_SCHWZ;
   double f_cut = -fgrav_pw( pl->M, r_cut );

   double f1 = 0.0;
   if( pw_flag ){
       f1 = -fgrav_pw( pl->M, script_r );
       if( f1 < f_cut ) //if f(r) more negative than f(1.5r_s)
            f1 = f_cut;
   }
   else{
       f1 = -fgrav( pl->M, script_r, pl->eps );
   }

   double cosa = dx/script_r;
   double sina = dy/script_r;

   double cosap = cosa*cosp+sina*sinp;
   double sinap = sina*cosp-cosa*sinp;

   if( mode==1 ){
      cosap = cosa*cos(pp)+sina*sin(pp);
      sinap = sina*cos(pp)-cosa*sin(pp);
   }
/*
   double rH = rp*pow( pl->M/3.,1./3.);
   double pd = 0.8; 
   double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));
*/

   double sint = script_r_perp/script_r;
   double cost = z/script_r;

   *fr = cosap*f1*sint; //*fd;
   *fp = sinap*f1*sint; //*fd;
   *fz = f1*cost;

}

void planet_src( struct planet * pl , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   double rp = xp[0];
   double rm = xm[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   
   double r = 0.5*(rp+rm);
   double vp  = r*omega;
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;
   double z = .5*(xp[2]+xm[2]);

   double Fr,Fp,Fz;
   planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 0 );

   cons[SRR] += rho*Fr*dVdt;
   cons[SZZ] += rho*Fz*dVdt;
   cons[LLL] += rho*Fp*r*dVdt;
   cons[TAU] += rho*( Fr*vr + Fz*vz + Fp*vp )*dVdt;

}

int nearest_planet_dist( struct domain *theDomain, double r, double phi ){
   
   int Npl = theDomain->Npl;
   double rmax = theDomain->theParList.rmax;
   struct planet *thePlanets = theDomain->thePlanets;

   int p;
   double nearest = 2*rmax;
   for( p=0; p<Npl; ++p ){
      struct planet *pl = thePlanets+p;
      double dr   = r - pl->r;
      double dp   = cos( phi - pl->phi );
      double dl2  = dr*dr + 2*r*(pl->r)*dp*dp;
      double dist = sqrt( dl2 );
      if( dist < nearest )
         nearest = dist;
   }
   return nearest;

}

void planet_RK_copy( struct planet * pl ){
   pl->RK_r     = pl->r;
   pl->RK_phi   = pl->phi;
   pl->RK_M     = pl->M;
   pl->RK_omega = pl->omega;
   pl->RK_vr    = pl->vr;
}

void planet_RK_adjust( struct planet * pl , double RK ){
   pl->r     = (1.-RK)*pl->r     + RK*pl->RK_r;
   pl->phi   = (1.-RK)*pl->phi   + RK*pl->RK_phi;
   pl->M     = (1.-RK)*pl->M     + RK*pl->RK_M;
   pl->omega = (1.-RK)*pl->omega + RK*pl->RK_omega;
   pl->vr    = (1.-RK)*pl->vr    + RK*pl->RK_vr;
}
