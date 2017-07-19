
#include "paul.h"
#include <string.h>

#define R_HOR 1.5

void initial( double * , double * );
double get_dV( double * , double * );
double get_dA( double *, double *, int);
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );
double get_moment_arm(double *, double *);
void subtract_omega( double * );
void reflect_prims(double *, double *, int);

void set_cell_init(struct cell *c, double *r_jph, double *z_kph, int j, int k)
{
    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
    double r = get_moment_arm(xp, xm);
    double phi = c->piph - 0.5*c->dphi;
    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
    initial(c->prim, x);
    subtract_omega(c->prim);
    if(NUM_C > BZZ)
    {
        double prim[NUM_Q], dA;
        if(NUM_FACES >= 1)
        {
            xm[1] = c->piph;
            x[1] = c->piph;
            dA = get_dA(xp, xm, 0);
            initial(prim, x);
            c->Phi[0] = prim[BPP]*dA;
        }
        if(NUM_FACES >= 3)
        {
            xm[1] = c->piph - c->dphi;
            x[1] = phi;

            xp[0] = r_jph[j-1];
            x[0] = r_jph[j-1];
            dA = get_dA(xp, xm, 1);
            initial(prim, x);
            c->Phi[1] = prim[BRR]*dA;
            xp[0] = r_jph[j];
            xm[0] = r_jph[j];
            x[0] = r_jph[j];
            dA = get_dA(xp, xm, 1);
            initial(prim, x);
            c->Phi[2] = prim[BRR]*dA;
        }
        if(NUM_FACES >= 5)
        {
            xm[0] = r_jph[j-1];
            x[0] = r;

            xp[2] = z_kph[j-1];
            x[2] = z_kph[j-1];
            dA = get_dA(xp, xm, 2);
            initial(prim, x);
            c->Phi[3] = prim[BZZ]*dA;
            xp[2] = z_kph[j];
            xm[2] = z_kph[j];
            x[2] = z_kph[j];
            dA = get_dA(xp, xm, 2);
            initial(prim, x);
            c->Phi[4] = prim[BZZ]*dA;
        }
    }
}

void set_cell_init_q(struct cell *c, double *r_jph, double *z_kph, 
                        int j, int k, int *qarr, int nq)
{
    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
    double r = get_moment_arm(xp, xm);
    double phi = c->piph - 0.5*c->dphi;
    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
    double temp_prim[NUM_Q];
    initial(temp_prim, x);
    subtract_omega(temp_prim);

    int iq;
    for(iq=0; iq<nq; iq++)
        c->prim[qarr[iq]] = temp_prim[qarr[iq]];
}

void set_cells_copy(struct cell *c, int Np, struct face *theFaces, 
                    int n0, int n1, int LR)
{
    int i,n,q;

    // Clear annulus.
    for(i=0; i<Np; i++)
    {
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] = 0.0;
        c[i].tempDoub = 0.0;
    }

    // Add outer strip.
    for(n=n0; n<n1; n++)
    {
        struct face *f = theFaces+n;
        struct cell *cDst, *cSrc;
        if(LR > 0)
        {
            cSrc = f->L;
            cDst = f->R;
        }
        else
        {
            cSrc = f->R;
            cDst = f->L;
        }
        for(q=0; q<NUM_Q; q++)
            (cDst->prim)[q] += (cSrc->prim)[q] * (f->dA);
        cDst->tempDoub += f->dA;
    }
    
    // Divide by total area.
    for(i=0; i<Np; i++)
    {
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] /= c[i].tempDoub;
    }
}

void set_cells_copy_distant(struct cell *c, int Np, struct cell *c1, int Np1)
{
    // Copies annulus c1 to c.

    int i, q;

    // Clear annulus.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] = 0.0;

    // Find first cell in c1 that intersects with c[0].
    int i1 = 0;
    double phi0 = c[0].piph - c[0].dphi;
    for(i1=0; i1 < Np1; i1++)
    {
        double dphip = c1[i1].piph - phi0;
        while(dphip > M_PI)
            dphip -= 2*M_PI;
        while(dphip < -M_PI)
            dphip += 2*M_PI;
        double dphim = dphip - c1[i1].dphi;

        if(dphip > 0 && dphim <= 0.0)
            break;
    }

    if(i1 >= Np1)
        i1 = 0;

    //Loop through all cells c, add contributions from all neighbours in c1.
    for(i=0; i<Np; i++)
    {
        double phip = c[i].piph;
        double phim = c[i].piph - c[i].dphi;
        double phip1, phim1;
        do
        {
            phip1 = c1[i1].piph;
            phim1 = c1[i1].piph - c1[i1].dphi;

            while(phip1 < phim)
            {
                phip1 += 2*M_PI;
                phim1 += 2*M_PI;
            }
            while(phim1 > phip)
            {
                phip1 -= 2*M_PI;
                phim1 -= 2*M_PI;
            }
            double phi1 = phip < phip1? phip : phip1;
            double phi2 = phim > phim1? phim : phim1;
            double dphi = phi1-phi2;

            if(dphi< 0.0)
                printf("WHOA %d %d %.12lg\n", i, i1, dphi);
            for(q=0; q<NUM_Q; q++)
                c[i].prim[q] += c1[i1].prim[q] * dphi;

            i1 = i1 < Np1-1 ? i1+1 : 0;
            if(i1 == Np1)
                i1 = 0;
        }
        while(phip1 < phip);
        i1 = i1 > 0 ? i1-1 : Np1-1;

    }

    // Divide by total area.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            c[i].prim[q] /= c[i].dphi;
}

void boundary_fixed_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<Ng; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-Ng; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<Ng; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}
void boundary_fixed_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-Ng; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
            }  
    }
}

void boundary_zerograd_rinn( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_1;
    int *fIndex = theDomain->fIndex_r;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=Ng-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                int JK = j+(Nr-1)*k;
                int n0 = fIndex[JK];
                int n1 = fIndex[JK+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, -1);

                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[URR] > 0)
                            theCells[jk][i].prim[URR] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_rout( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_1;
    int *fIndex = theDomain->fIndex_r;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-Ng; j<Nr; j++)
            {
                int jk = j + Nr*k;
                int JK = j-1 + (Nr-1)*k;
                int n0 = fIndex[JK];
                int n1 = fIndex[JK+1];
                
                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, +1);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[URR] < 0)
                            theCells[jk][i].prim[URR] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_zbot( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_2;
    int *fIndex = theDomain->fIndex_z;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j++)
            for(k=Ng-1; k>=0; k--)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk];
                int n1 = fIndex[jk+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, -1);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[UZZ] > 0)
                            theCells[jk][i].prim[UZZ] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_ztop( struct domain *theDomain, int diode)
{
    struct cell **theCells = theDomain->theCells;
    struct face *theFaces = theDomain->theFaces_2;
    int *fIndex = theDomain->fIndex_z;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-Ng; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int n0 = fIndex[jk-Nr];
                int n1 = fIndex[jk-Nr+1];

                set_cells_copy(theCells[jk], Np[jk], theFaces, n0, n1, +1);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        if(theCells[jk][i].prim[UZZ] < 0)
                            theCells[jk][i].prim[UZZ] = 0.0;
                    }
                }
            }
    }
}

void boundary_reflect_rinn( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=Ng-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*Ng-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
    }
}

void boundary_reflect_rout( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-Ng; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-Ng) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 0);
                }
            }
    }
}

void boundary_reflect_zbot( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=Ng-1; k>=0; k--)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*Ng-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
    }
}

void boundary_reflect_ztop( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-Ng; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-Ng) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(theCells[jk], Np[jk], 
                                        theCells[jk1], Np[jk1]);

                for(i=0; i<Np[jk]; i++)
                {
                    struct cell *c = &(theCells[jk][i]);
                    double xm[3] = {r_jph[j-1], c->piph - c->dphi, z_kph[k-1]};
                    double xp[3] = {r_jph[j  ], c->piph          , z_kph[k  ]};
                    double r = get_moment_arm(xp, xm);
                    double phi = c->piph - 0.5*c->dphi;
                    double x[3] = {r, phi, 0.5*(z_kph[k]+z_kph[k-1])};
                    reflect_prims(theCells[jk][i].prim, x, 2);
                }
            }
    }
}

void boundary_fixed_horizon( struct domain *theDomain)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int i,j,k;

    for(k=0; k<Nz; k++)
    {
        double zm = z_kph[k-1];
        double zp = z_kph[k];
        double zo = fabs(zp)>fabs(zm) ? zp : zm;

        if(fabs(zo) > R_HOR)
            continue;

        for(j=0; j<Nr; j++)
        {
            double ro = r_jph[j];

            double R = sqrt(zo*zo + ro*ro);
            
            if(R < R_HOR)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    set_cell_init(&(theCells[jk][i]), r_jph, z_kph, j, k);
                    theCells[jk][i].real = 0;
                }
            }  
        }
    }
}

void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<Ng; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}

void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-Ng; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}

void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<Ng; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k,
                                    q, nq);
            }  
    }
}
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq)
{
    struct cell **theCells = theDomain->theCells;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int Ng = theDomain->Ng;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-Ng; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(theCells[jk][i]), r_jph, z_kph, j, k, 
                                    q, nq);
            }  
    }
}
