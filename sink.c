#include "paul.h"

static int sinkType = 0;
static double sinkPar1 = 0.0;
static double sinkPar2 = 0.0;
static double sinkPar3 = 0.0;
static double sinkPar4 = 0.0;
static int nozzleType;
static double nozzlePar1 = 0.0;
static double nozzlePar2 = 0.0;
static double nozzlePar3 = 0.0;
static double nozzlePar4 = 0.0;
static double gamma_law = 0.0;
static int twoD = 0;

void setSinkParams(struct domain *theDomain)
{
    sinkType = theDomain->theParList.sinkType;
    sinkPar1 = theDomain->theParList.sinkPar1;
    sinkPar2 = theDomain->theParList.sinkPar2;
    sinkPar3 = theDomain->theParList.sinkPar3;
    sinkPar4 = theDomain->theParList.sinkPar4;
    nozzleType = theDomain->theParList.nozzleType;
    nozzlePar1 = theDomain->theParList.nozzlePar1;
    nozzlePar2 = theDomain->theParList.nozzlePar2;
    nozzlePar3 = theDomain->theParList.nozzlePar3;
    nozzlePar4 = theDomain->theParList.nozzlePar4;

    gamma_law = theDomain->theParList.Adiabatic_Index;
    if(theDomain->Nz == 1)
        twoD = 1;
}

double get_om(double *x);

void sink_src(double *prim, double *cons, double *xp, double *xm, double dVdt)
{
    if(nozzleType == 1)
    {
        double R0 = nozzlePar1;
        double dR = nozzlePar2;
        double v = nozzlePar3;
        double mach = nozzlePar4;

        double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);
        if(phi > M_PI) phi -= 2*M_PI;

        if((r-R0)*(r-R0) + r*r*phi*phi + z*z > dR*dR)
            return;

        double Mdot = 1.0;
        double rhodot = Mdot * 3.0/(4.0*M_PI*dR*dR*dR);
        if(twoD)
            rhodot *= 2*sqrt(dR*dR - (r-R0)*(r-R0) - r*r*phi*phi);
        if(rhodot * dVdt > cons[RHO])
            rhodot = cons[RHO] / dVdt;

        double x[3] = {r, phi, z};
        double om = get_om(x);
        double cs2 = v*v / (mach*mach);
        double eps = cs2 / (gamma_law*(gamma_law-1));

        cons[RHO] += rhodot * dVdt;
        cons[LLL] += rhodot * r*v * dVdt;
        cons[TAU] += (0.5* rhodot * (v-r*om)*(v-r*om) + rhodot*eps) * dVdt;
        if(NUM_Q > NUM_C)
            cons[NUM_C] += rhodot*dVdt;
    }

    if(sinkType == 2)
    {
     	double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);

        double cosp = cos(phi);
        double sinp = sin(phi);
        double gx = r*cosp
        double gy = r*sinp

        double px, py, pz, dx, dy, dz, mag, eps;
        double argTot = 0.0;
        //Pretend for now that p1 is an array with planet 1 coords, p2 is similar
        int pi;
        for (pi=0; pi<Npl; pi++){
            cosp = cos(thePlanets[pi].phi);
            sinp = sin(thePlanets[pi].phi);
            px = thePlanets[pi].r*cosp;
            py = thePlanets[pi].r*sinp;

            dx = gx-px;
            dy = gy-py;
            dz = z-thePlanets[pi].z;
            mag = dx*dx + dy*dy + dz*dz;
            mag = mag*mag;
            eps = thePlanets[pi].eps;
            eps = eps*eps*eps*eps;

            argTot += exp(-mag*mag/eps);
        }
        //sinkPar1 is the sink rate, should be > 1/100 (Duffell+ 2019)
        //sinkPar2 is the ratio floor

        double rate = sinkPar1*fmax(thePlanets[0].om,thePlanets[1].om);
        double surfdiff = rate*argTot;
        double ratio = 1.0-surfdiff;
        ratio = fmax(ratio,sinkPar2);
        CONS[URR] *= ratio;
        CONS[UZZ] *= ratio;
        CONS[UPP] *= ratio;
        CONS[RHO] *= ratio;
        CONS[TAU] *= ratio;
    }


}
