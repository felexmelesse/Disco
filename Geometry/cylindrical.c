
#include "../paul.h"

static double phi_max = 0.0;

void setGeometryParams( struct domain * theDomain ){
   phi_max = theDomain->phi_max;
}

double get_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp < 0.0 ) dp += phi_max;
   while( dp > phi_max) dp -= phi_max;
   return(dp);
}

double get_signed_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp <-.5*phi_max ) dp += phi_max;
   while( dp > .5*phi_max ) dp -= phi_max;
   return(dp);
}

double get_moment_arm( double * xp , double * xm ){
   double rp = xp[0];
   double rm = xm[0];
   double r2 = .5*(rp*rp+rm*rm);
   return( sqrt(r2) );
}

double get_centroid(double xp, double xm, int dim)
{
    if(dim == 1)
        return 2.0*(xp*xp+xp*xm+xm*xm) / (3.0*(xp+xm));
    else
        return 0.5*(xp+xm);
}

double get_dL( double * xp , double * xm , int dim ){
    double r = .5*(xp[0]+xm[0]);
    double dphi = get_dp(xp[1], xm[1]);
    if(dim == 0)
        return r*dphi;
    else if(dim == 1)
        return xp[0]-xm[0];
    else
        return xp[2]-xm[2];
}

double get_dA( double * xp , double * xm , int dim ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1], xm[1]);
    double dz   = xp[2]-xm[2];
    if(dim == 0)
        return dr*dz;
    else if(dim == 1)
        return r*dphi*dz;
    else
        return r*dr*dphi;
}

double get_dV( double * xp , double * xm ){
    double r  = .5*(xp[0]+xm[0]);
    double dr   = xp[0]-xm[0];
    double dphi = get_dp(xp[1],xm[1]);
    double dz   = xp[2]-xm[2];

    return( r*dr*dphi*dz );
}

double get_scale_factor( double * x, int dim)
{
    if(dim == 0)
        return x[0];
    return 1.0;
}

double get_vol_element(double *x)
{
    return x[0];
}

void get_xyz(double *x, double *xyz)
{
    xyz[0] = x[0] * cos(x[1]);
    xyz[1] = x[0] * sin(x[1]);
    xyz[2] = x[2];
}

void get_rpz(double *x, double *rpz)
{
    rpz[0] = x[0];
    rpz[1] = x[1];
    rpz[2] = x[2];
}

void get_coords_from_xyz(double *xyz, double *x)
{
    x[0] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
    x[1] = atan2(xyz[1],xyz[0]);
    x[2] = xyz[2];
}

void get_coords_from_rpz(double *rpz, double *x)
{
    x[0] = rpz[0];
    x[1] = rpz[1];
    x[2] = rpz[2];
}

void get_vec_rpz(double *x, double *v, double *vrpz)
{
    vrpz[0] = v[0];
    vrpz[1] = v[1];
    vrpz[2] = v[2];
}

void get_vec_from_rpz(double *x, double *vrpz, double *v)
{
    v[0] = vrpz[0];
    v[1] = vrpz[1];
    v[2] = vrpz[2];
}

void get_vec_xyz(double *x, double *v, double *vxyz)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    vxyz[0] = cp*v[0] - sp*v[1];
    vxyz[1] = sp*v[0] + cp*v[1];
    vxyz[2] = v[2];
}

void get_vec_from_xyz(double *x, double *vxyz, double *v)
{
    double phi = x[1];
    double cp = cos(phi);
    double sp = sin(phi);

    v[0] =  cp*vxyz[0] + sp*v[1];
    v[1] = -sp*vxyz[0] + cp*v[1];
    v[2] = vxyz[2];
}
