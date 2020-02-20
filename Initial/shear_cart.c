
#include "../paul.h"

static double nu = 0.0;
static double gam = 0.0;
static double x0 = 0.0;
static double cs0 = 0.0;
static double sig0 = 0.0;
static double v0 = 0.0;
static int prof_choice = 0;
static int isothermal_flag = 0;
static double *t = NULL;
    
void get_xyz(double *, double *);
void get_vec_from_xyz(double *, double *, double *);
void get_vec_contravariant(double *, double *, double *);

double get_cs2(double *x);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    nu = theDomain->theParList.viscosity;
    isothermal_flag = theDomain->theParList.isothermal_flag;
    prof_choice = theDomain->theParList.initPar0;
    x0 = theDomain->theParList.initPar1;
    cs0 = theDomain->theParList.initPar2;
    sig0 = theDomain->theParList.initPar3;
    v0 = theDomain->theParList.initPar4;
    t = &(theDomain->t);
}

void initial(double *prim, double *x)
{
    double xyz[3];
    get_xyz(x, xyz);
    double X = xyz[0];

    double rho, vx, vy, P;

    rho = 1.0;

    if(isothermal_flag)
        P = get_cs2(x) * rho / gam;
    else
        P = cs0*cs0 * rho;

    double t0 = sig0*sig0 / (2 * nu);
    double sig2 = 2 * nu * (*t+t0);

    vx = 0.0;
    vy = v0 * exp(-0.5*(X-x0)*(X-x0) / sig2) * sig0/sqrt(sig2);

    double Vxyz[3] = {vx, vy, 0};
    double V[3];
    get_vec_from_xyz(x, Vxyz, V);
    get_vec_contravariant(x, V, V);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = V[0];
    prim[UPP] = V[1];
    prim[UZZ] = V[2];

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
