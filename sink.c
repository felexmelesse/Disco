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

static int coolType = 0;
static double coolPar1 = 0.0;
static double coolPar2 = 0.0;
static double coolPar3 = 0.0;
static double coolPar4 = 0.0;

static double gamma_law = 0.0;
static int twoD = 0;
static int Npl = 0;
static struct planet *thePlanets = NULL;

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
    coolType = theDomain->theParList.coolType;
    coolPar1 = theDomain->theParList.coolPar1;
    coolPar2 = theDomain->theParList.coolPar2;
    coolPar3 = theDomain->theParList.coolPar3;
    coolPar4 = theDomain->theParList.coolPar4;

    gamma_law = theDomain->theParList.Adiabatic_Index;
    if(theDomain->Nz == 1)
        twoD = 1;

    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
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

    if (sinkType == 1) 		//sink a la Farris et al. 2014
    {


    }

    if(sinkType == 2)		//sink a la Duffell et al. 2019
    {
     	double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);

        double cosp = cos(phi);
        double sinp = sin(phi);
        double gx = r*cosp;
        double gy = r*sinp;

        double px, py, dx, dy, mag, eps;
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
            mag = dx*dx + dy*dy + z*z;
            mag = mag*mag;
            eps = thePlanets[pi].eps;
            eps = eps*eps*eps*eps;

            argTot += exp(-mag*mag/eps);
        }
        //sinkPar1 is the sink rate, should be > 1/100 (Duffell+ 2019)
        //sinkPar2 is the ratio floor

        double rate = sinkPar1*fmax(thePlanets[0].omega,thePlanets[1].omega);
        double surfdiff = rate*argTot;
        double ratio = 1.0-surfdiff;
        ratio = fmax(ratio,sinkPar2);
        cons[URR] *= ratio;
        cons[UZZ] *= ratio;
        cons[UPP] *= ratio;
        cons[RHO] *= ratio;
        cons[TAU] *= ratio;
    }
}

void cooling(double *prim, double *cons, double *xp, double *xm, double dt )
{
    if(coolType == 1)
    {
        double press = prim[PPP];
     	double r = 0.5*(xp[0]+xm[0]);
        double gm1 = gamma_law-1.0;
        double beta = coolPar1;
        double om = 1.0;
        if (r > 1.0) om = pow(r,-1.5);
        cons[TAU] = cons[TAU] - (press/gm1)*beta*om;
    }
    if(coolType == 2)
    {
        double tfact = coolPar1;	//sigma_sb/kappa?
        double press, T, gm1;
        press = prim[PPP];
        gm1 = gamma_law-1.0;
        T = press/prim[RHO];
        cons[TAU] -= (press/gm1)*(pow(1.0 + 8.0*gm1*tfact*T*T*T*dt/(prim[RHO]*prim[RHO]), -1.0/3.0) - 1.0);
    }
}
