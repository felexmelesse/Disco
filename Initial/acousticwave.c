
#include "../paul.h"

static double gam = 0.0;
static double rho_ref = 1.0;
static double cs_ref = 1.0;
static double mach = 0.0;
static double L = 0.0;
static double a = 0.0;
static double x0 = 0.0;
static double theta = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    mach = theDomain->theParList.initPar1;
    a = theDomain->theParList.initPar2;     //1.0 in RAM
    x0 = theDomain->theParList.initPar3;  
    L = theDomain->theParList.initPar4;     //0.3 in RAM
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];

    double X = (r*cos(phi-theta) - x0) / L;

    double rho, v, ur, up, P;
    
    double P_ref = rho_ref * cs_ref*cs_ref/gam;
    double v_ref = mach*cs_ref;

    double f = X*X-1.0;
    if(fabs(X) < 1.0)
        rho = rho_ref*(1.0 + a*f*f*f*f);
    else
        rho = rho_ref;

    P = P_ref * pow(rho/rho_ref, gam);

    double cs = sqrt(gam*P/rho);
    v = v_ref + 2*(cs-cs_ref)/(gam-1);

    ur = v * cos(phi-theta);
    up = -v * sin(phi-theta) / r;

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = ur;
    prim[UPP] = up;
    prim[UZZ] = 0.0;

    int q;
    for(q = 5; q < NUM_Q; q++)
        prim[q] = 0.0;
}
