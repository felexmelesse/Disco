#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double xi = 0.0;
static double delta = 0.0;
static double rho0 = 0.0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mdot = theDomain->theParList.initPar1;
    R0 = theDomain->theParList.initPar2;
    delta = theDomain->theParList.initPar5;
    xi = theDomain->theParList.initPar5;
    rho0 = theDomain->theParList.initPar3;
    Mach = theDomain->theParList.Disk_Mach;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain-> Npl;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double R = r/R0;
    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho, fact;

    double sig0 = Mdot*Mach*Mach/(3.0*3.14159265*visc);

    fact = exp(-pow((R),-xi));
    rho = sig0*pow(R,-delta)*fact+rho0;
    double dsdr = sig0*pow(r, -delta-1)*fact*(-delta);
    double om = 1.0;
    double rdodr = 1.0;
    double ddro = 0.0;
    if (r > 1.0){
      om = pow(r,3.0);
      om = om*pow((1.0 + 3.0/(r*r*16)), 2.0);
      dsdr = sig0*fact*pow(R, -delta)*(xi*pow(R, -xi - 1.0) - delta/R);
      double dpdr = dsdr*cs2 + rho/(r*r*Mach);
      om = sqrt(om + dpdr/(r*rho));
      rdodr = -0.5*pow(r, -1.5) - (15./32.)*R0*R0*pow(r, -3.5);
      ddro = 0.75*pow(r, -2.5) + (15.*7./64.)*pow(r, -4.5)*R0*R0;
    }
    double nu;
    if (Npl < 2){
        nu = visc*cs2/(om);
    }
    else{
       double omtot = 0;
       double cosp, sinp, px, py, dx, dy, gx, gy, mag;
       gx = r*sin(phi);
       gy = r*cos(phi);
       int np;
       for(np = 0; np<Npl; np++){
          cosp = cos(thePlanets[np].phi);
          sinp = sin(thePlanets[np].phi);
          px = thePlanets[np].r*cosp;
          py = thePlanets[np].r*sinp;
          dx = gx-px;
          dy = gy-py;
          mag = dx*dx + dy*dy + thePlanets[np].eps;
          omtot +=  thePlanets[np].M*pow(mag, -1.5);
       }
       nu = visc*cs2/sqrt(omtot);
    }

    double dnudr = -1.5*nu/r;
    double dvisc = r*r*rho*rdodr*dnudr + nu*r*r*rdodr*dsdr + rho*r*r*nu*ddro + 2*r*nu*rho*rdodr;

    double v = 2.0*dvisc/(r*r*om*rho);
    double P = rho*cs2/gam;
 
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)
        prim[q2] = 0.0;
}
