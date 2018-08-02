
#include "../paul.h"

static double gam = 0.0;
static double rho = 1.0;
static double cs = 1.0;
static double mach = 0.0;
static double costheta0 = 0.0;
static double phi0_o_pi = 0.0;

void get_vec_from_xyz(double *x, double *vxyz, double *v);
void get_vec_contravariant(double *x, double *v, double *vc);
    
void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    mach = theDomain->theParList.initPar1;
    cs = theDomain->theParList.initPar2;
    phi0_o_pi = theDomain->theParList.initPar3;
    costheta0 = theDomain->theParList.initPar4;
}

void initial( double * prim , double * x )
{
    double phi0 = phi0_o_pi * M_PI;
    double sintheta0 = sqrt(1.0-costheta0*costheta0);

    double kx = cos(phi0)*sintheta0;
    double ky = sin(phi0)*sintheta0;
    double kz = costheta0;

    double P = rho * cs*cs/gam;
    double v = mach*cs;

    double vxyz[3] = {v*kx, v*ky, v*kz};
    double u[3];
    get_vec_from_xyz(x, vxyz, u);
    get_vec_contravariant(x, u, u);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = u[0];
    prim[UPP] = u[1];
    prim[UZZ] = u[2];

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
