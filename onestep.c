
#include "paul.h"

void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clean_pi( struct domain * );
void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );
void adjust_RK_planets( struct domain * , double );
void move_cells( struct domain * , double );
void calc_dp( struct domain * );
void calc_prim( struct domain * );
void calc_cons( struct domain * );
void B_faces_to_cells( struct domain * , int );

void setup_faces( struct domain * , int );
void phi_flux( struct domain * , double dt );
void trans_flux( struct domain * , double dt , int );
void add_source( struct domain * , double dt );

void avg_Efields( struct domain * );
void update_B_fluxes( struct domain * , double );
void subtract_advective_B_fluxes( struct domain * );
void check_flipped( struct domain * , int );
void flip_fluxes( struct domain * , int );

void movePlanets( struct planet * , double , double );
int planet_motion_analytic(void);

void boundary_r( struct domain * );
void boundary_trans( struct domain * , int );
void exchangeData( struct domain * , int );

//int get_num_rzFaces( int , int , int );
int set_B_flag( void );

void checkNaNs(struct domain *theDomain, char label[])
{
    int Nz = theDomain->Nz;
    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;

    int i, j, k, q;
    int count_p = 0;
    int count_c = 0;
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell *c = &(theDomain->theCells[jk][i]);

                int flag = 0;
                for(q=0; q<NUM_Q; q++)
                {
                    if(c->prim[q] != c->prim[q])
                    {
                        count_p++;
                        flag = 1;

                    }
                    if(c->cons[q] != c->cons[q])
                    {
                        count_c++;
                        flag = 1;
                    }
                }
                //if(flag)
                //    printf("  NaN action at k = %d, j = %d, i = %d\n", k,j,i);
            }
        }
    if(count_p > 0)
        printf("NaNs in prim @ %s!\n", label);
    if(count_c > 0)
        printf("NaNs in cons @ %s!\n", label);
}

void print_Phi_border(struct domain *theDomain, int step)
{
    if(1)
        return;
    int *dim_rank = theDomain->dim_rank;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;

    int *fIr = theDomain->fIndex_r;
    int *fIz = theDomain->fIndex_z;
    struct face *fr = theDomain->theFaces_1;
    struct face *fz = theDomain->theFaces_2;
    double *zkph = theDomain->z_kph;
    double *rjph = theDomain->r_jph;
    double t = theDomain->t;

    FILE *F;
    char fname[256];

    int j,k;
    // R-
    sprintf(fname, "Phi_%.2e_%d_%02d_%02d_Ra.txt", t, step, dim_rank[0], dim_rank[1]);
    F = fopen(fname, "w");
    j = NgRa == 0 ? 0: NgRa-1;
    fprintf(F, "r = %.3lf\n", rjph[j]);
    for(k=0; k<Nz; k++)
    {
        int jk = j+(Nr-1)*k;
        int f;
        fprintf(F, "%03d %.3f %.3f:", k, zkph[k-1], zkph[k]);
        double PhiTot = 0.0;
        for(f = fIr[jk]; f < fIr[jk+1]; f++)
        {
            double Phi;
            if(fr[f].LRtype == 0)
                Phi = fr[f].L->Phi[2];
            else
                Phi = fr[f].R->Phi[1];
            PhiTot += Phi;
            fprintf(F, " %.12lg", Phi);
        }
        fprintf(F, " (%.12lg)\n", PhiTot);
    }
    fclose(F);

    // R+
    sprintf(fname, "Phi_%.2e_%d_%02d_%02d_Rb.txt", t, step, dim_rank[0], dim_rank[1]);
    F = fopen(fname, "w");
    j = NgRb == 0 ? Nr-2 : Nr-NgRb-1;
    fprintf(F, "r = %.3lf\n", rjph[j]);
    for(k=0; k<Nz; k++)
    {
        int jk = j+(Nr-1)*k;
        int f;
        fprintf(F, "%03d %.3f %.3f:", k, zkph[k-1], zkph[k]);
        double PhiTot = 0.0;
        for(f = fIr[jk]; f < fIr[jk+1]; f++)
        {
            double Phi;
            if(fr[f].LRtype == 0)
                Phi = fr[f].L->Phi[2];
            else
                Phi = fr[f].R->Phi[1];
            PhiTot += Phi;
            fprintf(F, " %.12lg", Phi);
        }
        fprintf(F, " (%.12lg)\n", PhiTot);
    }
    fclose(F);

    // Z-
    sprintf(fname, "Phi_%.2e_%d_%02d_%02d_Za.txt", t, step, dim_rank[0], dim_rank[1]);
    F = fopen(fname, "w");
    k = NgZa == 0 ? 0 : NgZa-1;
    fprintf(F, "z = %.3lf\n", zkph[k]);
    for(j=0; j<Nr; j++)
    {
        int jk = j+Nr*k;
        int f;
        fprintf(F, "%03d %.3f %.3f:", j, rjph[j-1], rjph[j]);
        double PhiTot = 0.0;
        for(f = fIz[jk]; f < fIz[jk+1]; f++)
        {
            double Phi;
            if(fr[f].LRtype == 0)
                Phi = fz[f].L->Phi[4];
            else
                Phi = fz[f].R->Phi[3];
            PhiTot += Phi;
            fprintf(F, " %.12lg", Phi);
        }
        fprintf(F, " (%.12lg)\n", PhiTot);
    }
    fclose(F);

    // Z+
    sprintf(fname, "Phi_%.2e_%d_%02d_%02d_Zb.txt", t, step, dim_rank[0], dim_rank[1]);
    F = fopen(fname, "w");
    k = NgZb == 0 ? Nz-2 : Nz-NgZb-1;
    fprintf(F, "z = %.3lf\n", zkph[k]);
    for(j=0; j<Nr; j++)
    {
        int jk = j+Nr*k;
        int f;
        fprintf(F, "%03d %.3f %.3f:", j, rjph[j-1], rjph[j]);
        double PhiTot = 0.0;
        for(f = fIz[jk]; f < fIz[jk+1]; f++)
        {
            double Phi;
            if(fr[f].LRtype == 0)
                Phi = fz[f].L->Phi[4];
            else
                Phi = fz[f].R->Phi[3];
            PhiTot += Phi;
            fprintf(F, " %.12lg", Phi);
        }
        fprintf(F, " (%.12lg)\n", PhiTot);
    }
    fclose(F);

    sprintf(fname, "Grid_%.2e_%d_%02d_%02d.txt", t, step, dim_rank[0],
                dim_rank[1]);
    F = fopen(fname, "w");
    for(k=0; k<Nz; k++)
    {
        fprintf(F, "%03d %.3f %.3f:\n", k, zkph[k-1], zkph[k]);
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            fprintf(F, "    %03d %.3f %.3f:", j, rjph[j-1], rjph[j]);
            int i;
            for(i=0; i<theDomain->Np[jk]; i++)
                fprintf(F, " %.3lf", theDomain->theCells[jk][i].piph);
            fprintf(F, "\n");
        }
    }
    fclose(F);

    sprintf(fname, "Grid_Phi_%.2e_%d_%02d_%02d.txt", t, step, dim_rank[0],
                dim_rank[1]);
    F = fopen(fname, "w");
    for(k=0; k<Nz; k++)
    {
        fprintf(F, "%03d %.3f %.3f:\n", k, zkph[k-1], zkph[k]);
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            fprintf(F, "    %03d %.3f %.3f:\n", j, rjph[j-1], rjph[j]);
            int i;
            for(i=0; i<theDomain->Np[jk]; i++)
            {
                fprintf(F, "        %d %.3lf:", i, 
                        theDomain->theCells[jk][i].piph);
                int q;
                for(q=0; q<NUM_FACES; q++)
                    fprintf(F, " %.12lg", theDomain->theCells[jk][i].Phi[q]);
                fprintf(F, "\n");
            }
        }
    }
    fclose(F);

    sprintf(fname, "Grid_B_%.2e_%d_%02d_%02d.txt", t, step, dim_rank[0],
                dim_rank[1]);
    F = fopen(fname, "w");
    for(k=0; k<Nz; k++)
    {
        fprintf(F, "%03d %.3f %.3f:\n", k, zkph[k-1], zkph[k]);
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            fprintf(F, "    %03d %.3f %.3f:\n", j, rjph[j-1], rjph[j]);
            int i;
            for(i=0; i<theDomain->Np[jk]; i++)
            {
                fprintf(F, "        %d %.3lf:", i, 
                        theDomain->theCells[jk][i].piph);
                fprintf(F, " %.12lg", theDomain->theCells[jk][i].cons[BRR]);
                fprintf(F, " %.12lg", theDomain->theCells[jk][i].cons[BPP]);
                fprintf(F, " %.12lg", theDomain->theCells[jk][i].cons[BZZ]);
                fprintf(F, "\n");
            }
        }
    }
    fclose(F);

    double phir[Nz][Nr-1];
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr-1; j++)
            phir[k][j] = 0.0;
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            int i;
            for(i=0; i<theDomain->Np[jk]; i++)
            {
                if(j > 1)
                    phir[k][j-1] += theDomain->theCells[jk][i].Phi[1];
                if(j < Nr-1)
                    phir[k][j] += theDomain->theCells[jk][i].Phi[2];
            }
        }

    double phiz[Nz-1][Nr];
    for(k=0; k<Nz-1; k++)
        for(j=0; j<Nr; j++)
            phiz[k][j] = 0.0;
    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = j + Nr*k;
            int i;
            for(i=0; i<theDomain->Np[jk]; i++)
            {
                if(k > 1)
                    phir[k-1][j] += theDomain->theCells[jk][i].Phi[3];
                if(k < Nz-1)
                    phir[k][j] += theDomain->theCells[jk][i].Phi[4];
            }
        }

    sprintf(fname, "Phi_R_%.2e_%d_%02d_%02d.txt", t, step, dim_rank[0],
                dim_rank[1]);
    F = fopen(fname, "w");
    for(k=0; k<Nz; k++)
    {
        fprintf(F, "%03d %.3f %.3f:\n", k, zkph[k-1], zkph[k]);
        for(j=0; j<Nr-1; j++)
            fprintf(F, " %.12lg", phir[k][j]);
        fprintf(F, "\n");
    }
    fclose(F);

    sprintf(fname, "Phi_Z_%.2e_%d_%02d_%02d.txt", t, step, dim_rank[0],
                dim_rank[1]);
    F = fopen(fname, "w");
    for(k=0; k<Nz-1; k++)
    {
        fprintf(F, "%03d %.3f %.3f:\n", k, zkph[k-1], zkph[k]);
        for(j=0; j<Nr; j++)
            fprintf(F, " %.12lg", phiz[k][j]);
        fprintf(F, "\n");
    }
    fclose(F);
}


void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step , double global_dt ){

   int Nz = theDomain->Nz;
   int bflag = set_B_flag();

   if( first_step ) set_wcell( theDomain );
   adjust_RK_cons( theDomain , RK );

   phi_flux( theDomain , dt );

   setup_faces( theDomain , 1 );
   trans_flux( theDomain , dt , 1 );

   if( Nz > 1 ){
      setup_faces( theDomain , 2 );
      trans_flux( theDomain , dt , 2 );
   }

   print_Phi_border(theDomain, 0);

   if( bflag && NUM_EDGES >= 4 ){
      avg_Efields( theDomain );
      subtract_advective_B_fluxes( theDomain );
      update_B_fluxes( theDomain , dt );
   }
   
   print_Phi_border(theDomain, 1);

   add_source( theDomain , dt );

   if( first_step ){
      move_cells( theDomain , dt );
      if( bflag ){
         check_flipped( theDomain , 0 );
         flip_fluxes( theDomain , 0 );
         if( Nz>1 ){
            check_flipped( theDomain , 1 );
            flip_fluxes( theDomain , 1 );
         }
      }
   }
   

   if( !planet_motion_analytic() || first_step ){
      adjust_RK_planets( theDomain , RK );
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
   }
   clean_pi( theDomain );
   calc_dp( theDomain );
   
   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 1 );
   }
   
   print_Phi_border(theDomain, 2);

   calc_prim( theDomain ); //ORDERING??? AFTER?
   
   /*if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 0 );
   }

   //TODO: interaction with MHD? Hail Mary
   calc_cons(theDomain);
   */

   exchangeData( theDomain , 0 );
   if(! theDomain->theParList.R_Periodic)
      boundary_trans( theDomain , 1 );
   if( Nz > 1 ){
      exchangeData( theDomain , 1 );
      if(! theDomain->theParList.Z_Periodic)
         boundary_trans( theDomain , 2 );
   }
   
   print_Phi_border(theDomain, 3);

   //TODO: This was BEFORE BCs, but if wrecks cell pointers...
   //      Here, the BCs may not be satisfied if boundary zones are AMR'd...
   //TODO 2: AMR leading to STRANGE behaviour in 3d? Z boundaries being
   //           overwritten?  Needs a closer look.
   if( last_step ){
      //AMR( theDomain );
   }


   if( theDomain->theFaces_1 ) free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );

}
