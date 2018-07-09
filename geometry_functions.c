
#include "paul.h"

double get_centroid(double xp, double xm, int dim);
double get_scale_factor( double * x, int dim);

void get_centroid_arr(double *xp, double *xm, double *x)
{
    x[0] = get_centroid(xp[0], xm[0], 1);
    x[1] = get_centroid(xp[1], xm[1], 0);
    x[2] = get_centroid(xp[2], xm[2], 2);
}

void get_vec_contravariant(double *x, double *v, double *vc)
{
    vc[0] = v[0] / get_scale_factor(x, 1);
    vc[1] = v[1] / get_scale_factor(x, 0);
    vc[2] = v[2] / get_scale_factor(x, 2);
}

void get_vec_covariant(double *x, double *v, double *vc)
{
    vc[0] = v[0] * get_scale_factor(x, 1);
    vc[1] = v[1] * get_scale_factor(x, 0);
    vc[2] = v[2] * get_scale_factor(x, 2);
}
