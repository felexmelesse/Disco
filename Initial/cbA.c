
#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double rho0 = 0.0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mdot = theDomain->theParList.initPar1;
    R0 = theDomain->theParList.initPar2;
    rho0 = theDomain->theParList.initPar3;
    Mach = theDomain->theParList.Disk_Mach;
    thePlanets = theDomain->thePlanets;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double R = r/R0;
    double phi = x[1];

    double cs2 = get_cs2(x);

    double xi = 4.0;
    double p = 0.5;
    double rho;

    double sig0 = Mdot*Mach*Mach/(3.0*3.14159265*visc);

    rho = sig0*pow(R,-p)*exp(-pow((R+0.5),-xi))+0.01;
    double om = 1.0;
    if (r>1.0) om = pow(r,-1.5);

    double nu;
    nu = visc*cs2/om;

    double v = -1.5*nu/r;
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
