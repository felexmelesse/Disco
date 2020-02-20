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

static double visc = 0.0;
static double Mach = 1.0;

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

    visc = theDomain->theParList.viscosity;
    Mach = theDomain->theParList.Disk_Mach;

    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
}

double get_om(double *x);
double get_cs2(double *);

void sink_src(double *prim, double *cons, double *xp, double *xm, double dV, double dt)
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
        if(rhodot * dV*dt > cons[RHO])
            rhodot = cons[RHO] / dV*dt;

        double x[3] = {r, phi, z};
        double om = get_om(x);
        double cs2 = v*v / (mach*mach);
        double eps = cs2 / (gamma_law*(gamma_law-1));

        cons[RHO] += rhodot * dV*dt;
        cons[LLL] += rhodot * r*v * dV*dt;
        cons[TAU] += (0.5* rhodot * (v-r*om)*(v-r*om) + rhodot*eps) * dV*dt;
        if(NUM_Q > NUM_C)
            cons[NUM_C] += rhodot*dV*dt;
    }

    //sink a la Farris et al. 2014
    if (sinkType == 1)
    {
      double r = 0.5*(xp[0]+xm[0]);
      double phi = 0.5*(xp[1]+xm[1]);
      double z = 0.5*(xp[2]+xm[2]);

      double cosp = cos(phi);
      double sinp = sin(phi);
      double gx = r*cosp;
      double gy = r*sinp;

      double px, py, dx, dy, mag, eps;
      double arg = 0.0;
      double argt = 0.0;
      double factor = Mach*Mach/(3.0*M_PI*visc*thePlanets[0].omega);		//assumes alpha viscosity for now
      int pi;
      for (pi=0; pi<Npl; pi++)
      {
          cosp = cos(thePlanets[pi].phi);
          sinp = sin(thePlanets[pi].phi);
          px = thePlanets[pi].r*cosp;
          py = thePlanets[pi].r*sinp;
          dx = gx-px;
          dy = gy-py;
          eps = thePlanets[pi].eps;
          mag = sqrt(dx*dx + dy*dy + z*z + eps*eps);

          if (mag < 0.5)
          {
            arg = factor*pow(mag, 1.5)*pow(thePlanets[pi].M, -0.5);

            if (arg>0)
            {
               argt += dV/(dt*arg); 
               thePlanets[pi].dM += prim[RHO]*dV/(dt*arg);
            }
          }          
      }
      double ratio = 1.0;
      ratio = fmax(1.0 - argt, 1.e-5);

      cons[RHO] *= ratio;
      cons[URR] *= ratio;
      cons[UZZ] *= ratio;
      cons[UPP] *= ratio;
      cons[TAU] *= ratio;
    }

    //sink a la Duffell et al. 2019
    if(sinkType == 2)		
    {
     	double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);

        double rho = prim[RHO];
        double vr  = prim[URR];
        double vp  = prim[UPP]*r;
        double vz  = prim[UZZ];

        //should really evaluate at center
        double cs2m = get_cs2(xm);
        double cs2p = get_cs2(xp);
        double cs2 = 0.5*(cs2m + cs2p);

        double cosp = cos(phi);
        double sinp = sin(phi);
        double gx = r*cosp;
        double gy = r*sinp;

        double px, py, dx, dy, mag, eps;
        double rate, surfdiff;
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

            double arg = exp(-mag/eps);
            rate = sinkPar1*thePlanets[pi].omega;
            surfdiff = rho*rate*arg;
            thePlanets[pi].dM += surfdiff*dV;
            cons[DDD] -= surfdiff*dV*dt;
            cons[TAU] -= surfdiff*cs2*dV*dt;
            cons[SZZ] -= vz*surfdiff*dV*dt;
            if (sinkPar2 == 0.0){
              cons[SRR] -= vr*surfdiff*dV*dt;
              cons[LLL] -= r*vp*surfdiff*dV*dt;
             }
            else {
              double delta = fmin(sinkPar2, 1.0);
              double rp, omp, vxp, vyp, vxg, vyg, vxr, vyr, vp_p, vp_r, thetap, thetag, vxn, vyn, cphi, sphi, vg_r, vg_p;
              rp = thePlanets[pi].r;
              omp = thePlanets[pi].omega;

              vp_p = rp*omp;
              vp_r = thePlanets[pi].vr;
              thetap = 0.5*M_PI - thePlanets[pi].phi;
              vxp = vp_p*cos(thetap) - vp_r*sin(thetap);
              vyp = vp_p*sin(thetap) + vp_r*cos(thetap);
              
              thetag = 0.5*M_PI - phi;
              vxg = vp*cos(thetag) - vr*sin(thetag);
              vyg = vp*sin(thetag) + vr*cos(thetag);
             
              vxr = vxg - vxp;
              vyr = vyg - vyp;
              //dphi = atan2(dy, dx);
              cphi = dx/pow(mag,0.25);
              sphi = dy/pow(mag,0.25);
              vxn = vxr*( (1-delta)*cphi*cphi + sphi*sphi) + delta*vyr*sphi*cphi;
              vyn = vxr*( delta*cphi*sphi ) + vyr*(cphi*cphi + (1-delta)*sphi*sphi);

              vxg = vxn + vxp;
              vyg = vyn + vyp;
              vg_r = vxg*cos(thetag) + vyg*sin(thetag);
              vg_p = -vxg*sin(thetag) + vyg*cos(thetag);

              cons[SRR] -= vg_r*surfdiff*dV*dt;
              cons[LLL] -= r*vg_p*surfdiff*dV*dt;
            }
        }
    }
}

void cooling(double *prim, double *cons, double *xp, double *xm, double dV, double dt )
{
    if(coolType == 1)
    {
        double press = prim[PPP];
     	double r = 0.5*(xp[0]+xm[0]);
        double gm1 = gamma_law-1.0;
        double beta = coolPar1;
        double om = 1.0;
        if (r > 1.0) om = pow(r,-1.5);
        cons[TAU] -= (press/gm1)*beta*om*dt*dV;
        //cons[TAU] -= (press/gm1)*dt*dV;
        
    }
    if(coolType == 2)
    {
        double beta = coolPar1;	
        double press, gm1;
        double sigma = prim[RHO];
        press = prim[PPP];
        gm1 = gamma_law-1.0;
        //T = press/sigma;
        double arg = 1.0 + 8*gm1*beta*press*press*press*dt/(sigma*sigma*sigma*sigma*sigma);
        cons[TAU] += (press/gm1)*dV*(pow(arg, -1./3.) - 1.0);
    }
    if(coolType == 3)
    {
        //constant H/r, assumes visc -> alpha
     	double r = 0.5*(xp[0]+xm[0]);
        double press, gm1;
        double sigma = prim[RHO];
        press = prim[PPP];
        gm1 = gamma_law-1.0;

        //here visc = alpha
        //double arg = 1.0 + 8*gm1*dt*(sigma/press)*(27./32.)*visc*beta*beta*r*r*pow(omega, 3.0);
        double arg = 1.0 + (27./8.)*gm1*dt*visc*Mach*Mach*pow(press, 3.0)*pow(sigma, -4.0); //Farris et al. 2015
        cons[TAU] += (press/gm1)*dV*(pow(arg, -1./3.) - 1.0);
    }
}
