
#include "../paul.h"
#include "../Calc/calc.h"

static double M = 0.0;
static double gam = 0.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
   M = 1.0;
}

void initial( double * prim , double * x ){

    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z*z);
    double R0 = 10.0;
    double rs = 10.0;
    double Mdot = 1.0;
    double b0 = 1.0e-4;
    
    double us2 = M / (2*rs);
    double as2 = us2;
    double us = -sqrt(us2);

    double rhos = -Mdot / (4*M_PI*rs*rs*us);
    double K = as2 / (gam * pow(rhos,gam-1));

    double a02 = 0.5*(5.0-3*gam)*as2;
    double rho0 = pow(a02 / (gam * K), 1.0/(gam-1));


    double uR = bondi_newt_solve(Mdot, M, R, gam, rho0, sqrt(a02), K);
    double rhoB = -Mdot / (4*M_PI*R*R*uR);
    double PB = K * pow(rhoB, gam);

    double chi = R>R0 ? 1.0 : 0.5*(1+tanh(tan(M_PI*(R/R0-0.5))));

    double rho = rhoB;
    double P = PB;

    if(R < rs)
    {
        rho = rhos;
        P = K * pow(rhos, gam);
        uR = -sqrt(2*M/R);
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = r/R * uR * chi;
    prim[UPP] = 0.0;
    prim[UZZ] = z/R * uR * chi;

    if(NUM_C > 5)
    {
        prim[BRR] = r/(R*R*R)*b0 * chi;
        prim[BPP] = 0.0;
        prim[BZZ] = z/(R*R*R)*b0 * chi;
    }

    int q;
    for(q=NUM_C; q < NUM_Q; q++)
    {
        prim[q] = 0.0;
    }
}
