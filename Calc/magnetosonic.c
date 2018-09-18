#include <math.h>
#include "integrate.h"

double ms_f(double x, void *args)
{
    double cs20 = ((double *)args)[0];
    double cA20 = ((double *)args)[1];
    double gam =  ((double *)args)[2];
    double cs2 = cs20*pow(x, gam-1.0);
    double cA2 = cA20*x;
    return sqrt(cs2 + cA2)/x;
}

double calc_magnetosonic_cf_int_newt(double rho, double rho0, double cs0, 
                                        double cA0, double gam)
{
    double args[3] = {cs0*cs0, cA0*cA0, gam};
    double cf = sqrt(cs0*cs0*pow(rho/rho0, gam-1.0) + cA0*cA0*rho/rho0);

    double intcf = romb(ms_f, 1.0, rho/rho0, 1025, 1.0e-12*cf, 1.0e-12, args);

    return intcf;
}
