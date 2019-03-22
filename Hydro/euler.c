
#include "../paul.h"

double get_om( double );
double get_om1( double );
double get_cs2( double, double, struct planet * );

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static double explicit_viscosity = 0.0;
static int include_viscosity = 0;
static int isothermal = 0;
static int alpha_flag = 0;

//Sink Params
static int    sink_flag = 0;
static double R_SINK    = 0.0;
static double TAU_SINK  = 1e-4;  //5.0;
static double pefficiency = 1.0;

static double a = 1.0;
static double q_planet = 1.0;
static int cs2Choice = 0;
static double G_EPS = -0.02;
static double t_step = 0.0;

//Variables for boundary damping
static double r_damp = 80.0;
static double r_max  = 0.0;

//Variables for calculating torques
double torque_accr = 0.0;
double torque_grav = 0.0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
   explicit_viscosity = theDomain->theParList.viscosity;
   include_viscosity = theDomain->theParList.visc_flag;
   alpha_flag = theDomain->theParList.alpha_flag;

//========= Params for Sink  =========================
   sink_flag = theDomain->theParList.sink_flag;
   R_SINK = theDomain->theParList.r_sink;
   TAU_SINK = theDomain->theParList.t_sink;
   pefficiency = theDomain->theParList.p_eff;

//========= Params for Binary Sims ===================
   cs2Choice = theDomain->theParList.Cs2_Profile;
   q_planet = theDomain->theParList.Mass_Ratio;
   G_EPS = theDomain->theParList.g_eps;

//========= Params for Damping =======================
    r_max = theDomain->theParList.rmax;
    //r_damp = theDomain->theParList.r_damp;

}

int set_B_flag(void){
   return(0);
}

void set_hydro_time( double t ){
    t_step = t;
}

double get_omega( double * prim , double * x ){
   return( prim[UPP] );
}

double Theta( double r ){
    
    double theta = (r - r_damp)/(r_max - r_damp);
    return theta*theta;
}

void planetaryForce( struct planet * , double , double, double , double * , double * , double*, int );

void prim2cons( double * prim , double * cons , double * x , double dV ){

   double r = x[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( r );
   double vp_off = vp - om*r;

   double v2  = vr*vr + vp_off*vp_off + vz*vz;

   double rhoe = Pp/(gamma_law - 1.);

   cons[DDD] = rho*dV;
   cons[TAU] = (.5*rho*v2 + rhoe )*dV;
   cons[SRR] = rho*vr*dV;
   cons[LLL] = r*rho*vp*dV;
   cons[SZZ] = rho*vz*dV;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      cons[q] = prim[q]*cons[DDD];
   }
}

void getUstar( double * prim , double * Ustar , double * x , double Sk , double Ss , double * n , double * Bpack ){

   double r = x[0];
   double rho = prim[RHO];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double Pp  = prim[PPP];

   double om = get_om( r );
   double vp_off = vp - om*r;
   double v2 = vr*vr+vp_off*vp_off+vz*vz;

   double vn = vr*n[0] + vp*n[1] + vz*n[2];

   double vn_off = vn - om*r*n[1];
   double Ss_off = Ss + vn_off - vn;

   double rhoe = Pp/(gamma_law - 1.);

   double rhostar = rho*(Sk - vn)/(Sk - Ss);
   double Pstar = Pp*(Ss - vn)/(Sk - Ss);
   double Us = rhoe*(Sk - vn)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] =   rhostar*( vr + (Ss-vn)*n[0] );
   Ustar[LLL] = r*rhostar*( vp + (Ss-vn)*n[1] );
   Ustar[SZZ] =   rhostar*( vz + (Ss-vn)*n[2] );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss_off*(Ss - vn) + Pstar;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void cons2prim( double * cons , double * prim , double * x , double dV, struct planet * thePlanets ){
   
   double r = x[0];
   double phi = x[1];
   double rho = cons[DDD]/dV;
   if( rho < RHO_FLOOR )   rho = RHO_FLOOR;
   double Sr  = cons[SRR]/dV;
   double Sp  = cons[LLL]/dV/r;
   double Sz  = cons[SZZ]/dV;
   double E   = cons[TAU]/dV;
   double om  = get_om( r );
   
   double vr = Sr/rho;
   double vp = Sp/rho;
   double vp_off = vp - om*r;
   double vz = Sz/rho;

   double KE = .5*( Sr*vr + rho*vp_off*vp_off + Sz*vz );
   double rhoe = E-KE;
   double Pp = (gamma_law - 1.)*rhoe;

   
//==================== Yike Cap Stuff ==================
   double CS_CAP = 10.0;
   double CS_FLOOR = 1e-5;
   double VEL_CAP = 20.0;
  
   if( isothermal ){
      double cs2 = get_cs2( r, phi, thePlanets );
      if( cs2 < CS_FLOOR*CS_FLOOR )
          cs2 = CS_FLOOR*CS_FLOOR;
      if( cs2 > CS_CAP*CS_CAP )
          cs2 = CS_CAP*CS_CAP;
      Pp = cs2*rho/gamma_law;
   }

   double v_mag = sqrt( 2.*KE/rho );
   if( v_mag > VEL_CAP ){
        vr = vr*VEL_CAP/v_mag;
        vp = vp*VEL_CAP/v_mag;
        vz = vz*VEL_CAP/v_mag;
   }
//========================================================

/*
   if( Pp  < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;
   if( isothermal ){
      double cs2 = get_cs2( r, phi );
      Pp = cs2*rho/gamma_law;
   }
*/

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = vp/r;
   prim[UZZ] = vz;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void flux( double * prim , double * flux , double * x , double * n ){
   
   double r = x[0];
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[URR];
   double vp  = prim[UPP]*r;
   double vz  = prim[UZZ];
   double om  = get_om( r );

   double vn = vr*n[0] + vp*n[1] + vz*n[2];
   double wn = om*r*n[1];
   double vp_off = vp - om*r;

   double rhoe = Pp/(gamma_law - 1.);
   double v2 = vr*vr + vp_off*vp_off + vz*vz;

   flux[DDD] = rho*vn;
   flux[SRR] = rho*vr*vn + Pp*n[0];
   flux[LLL] = r*rho*vp*vn + r*Pp*n[1];
   flux[SZZ] = rho*vz*vn + Pp*n[2];
   flux[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn - Pp*wn;

   int q;
   for( q=NUM_C ; q<NUM_Q ; ++q ){
      flux[q] = prim[q]*flux[DDD];
   }
   
}

double nearest_planet_dist( struct domain *, double, double );
double get_moment_arm( double * , double * );
double get_dp( double , double );
double get_planet_2dist( double, double, double, double );

void source( struct domain *theDomain, double * prim , double * cons , double * xp , double * xm , double dV, double dt, int last_step ){
  
   double rp = xp[0];
   double rm = xm[0];
   double dphi = get_dp(xp[1],xm[1]);
   double phip = xp[1];
   double phim = phip - dphi;
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double r_1  = .5*(rp+rm);
   double r2_3 = (rp*rp + rp*rm + rm*rm)/3.;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   double dVdt = dV*dt;


   double centrifugal = rho*omega*omega*r2_3/r_1*sin(.5*dphi)/(.5*dphi);
   double press_bal   = Pp/r_1;

   cons[SRR] += dVdt*( centrifugal + press_bal ); 

   double om  = get_om( r_1 );
   double om1 = get_om1( r_1 );

   cons[TAU] += dVdt*rho*vr*( om*om*r2_3/r_1 - om1*(omega-om)*r2_3 );

//================ Yike Binary Visc. Calc. =================
//TODO:   For eccentric or live binaries, need to fix this
//          -> source() has acces to planet objects

   double phi = 0.5*(phip+phim);
   double eps = G_EPS;
/*
   double om_bh = pow( a, -1.5 );
   double mu = q_planet/(1.+q_planet);

   double M1 = 1.-mu;
   double r1 = a*mu;
   double phi1 = om_bh*t_step - M_PI;
   double d1 = get_planet_2dist( r_1, phi, r1, phi1 );

   double M2 = mu;
   double r2 = a*(1.-mu);
   double phi2 = om_bh*t_step;
   double d2 = get_planet_2dist( r_1, phi, r2, phi2 );
*/
   struct planet *thePlanets = theDomain->thePlanets;
   double M1 = thePlanets[0].M;
   double M2 = thePlanets[1].M;
   double d1 = get_planet_2dist( r_1, phi, thePlanets[0].r, thePlanets[0].phi );
   double d2 = get_planet_2dist( r_1, phi, thePlanets[1].r, thePlanets[1].phi );

//==========================================================
 
   if( include_viscosity ){
      double nu = explicit_viscosity;
      if( alpha_flag ){
         double alpha = explicit_viscosity;
         double c = sqrt( gamma_law*prim[PPP]/prim[RHO] );
         double h = c*pow( r_1 , 1.5 );
        //viscosity calculation for Binary system--based on Yike's
         if( cs2Choice==4 ){
            double denom = 1./sqrt( M1*pow(d1*d1+eps*eps,-1.5) + M2*pow(d2*d2+eps*eps,-1.5) );
            h = c*denom;
         }
         nu = alpha*c*h;
      }
      //Seems like Yike doesn't have this source term???
      //Commented out just for test
      cons[SRR] += -dVdt*nu*rho*vr/(r_1*r_1);  
   }

}

void density_sink( struct domain *theDomain, double *prim, double *cons, double *xp, double *xm, double dV, double dt ){

//========================== This function is outdated ============================
//                              -> density_sink_yike() is the standard now

    double rho = prim[RHO];
    double vr  = prim[URR];
    double om  = prim[UPP];
    double vz  = prim[UZZ];
    double rp = xp[0];
    //double rm = xm[0];
    double r = get_moment_arm( xp, xm );
    //double r = 0.5*(xp[0]+xm[0]);
    double dphi = get_dp( xp[1],xm[1] );
    double phip = xp[1];
    double phim = phip - dphi;
    double dist = nearest_planet_dist( theDomain, rp, 0.5*(phip+phim) );
    if( sink_flag && dist < R_SINK ){
        //rho -= ( dt/TAU_SINK )*rho;
        //rho -= (1e-10)*rho;
        double sink_frac = dt/TAU_SINK;
        if( sink_frac > 1.0 ){
            sink_frac = 0.95;
            //printf("WARNING: TAU_SINK < dt \n");
        }
        //printf("sink_frac %f\n", sink_frac);
        prim[RHO] -= sink_frac*rho;
        //cons[DDD] -= sink_frac*(rho*dV);
        //cons[SRR] -= sink_frac*(rho*vr*dV);
        //cons[LLL] -= sink_frac*(rho*om*r*r*dV);
        //cons[SZZ] -= sink_frac*(rho*vz*dV);
        int q;
        double floor = 1e-5;
        for( q=0; q<NUM_Q; q++ ){
          cons[q] -= sink_frac*cons[q];
        }
        //double x[3] = {r, 0.5*(phip+phim), 0.5*(xp[2]+xm[2])};
        //prim2cons( prim, cons, x, dV );
    }
}

void density_sink_yike( struct domain *theDomain, double *prim, double *cons, double *xp, double *xm, double dV, double dt ){
    //This one is currently the sink being used (called by source() in misc.c)
    
    double rho = prim[RHO];
    double vr  = prim[URR];
    double om  = prim[UPP];
    double vz  = prim[UZZ];
    double rp = xp[0];
    double rm = xm[0];
    //double r = get_moment_arm( xp, xm );
    double r = 0.5*(xp[0]+xm[0]);
    double dphi = get_dp( xp[1],xm[1] );
    double phip = xp[1];
    double phim = phip - dphi;
    double phi  = 0.5*(phip+phim);
    //double dist = nearest_planet_dist( theDomain, 0.5*(rp+rm), 0.5*(phip+phim) );
    //TODO: Need to modify this to also save the nearest planet pointer
    //       For tracking accretion: theBH = thePlanet_nearest

    double *torques = theDomain->torques;

    int Npl = theDomain->Npl;
    struct planet *thePlanets = theDomain->thePlanets;
    
    int p;
    for( p=0; p<Npl; ++p ){
        struct planet *pl = thePlanets+p;
        double rbh = pl->r;
        double vbh = rbh*(pl->omega);
        double dist = get_planet_2dist( r, phi, pl->r, pl->phi);
        if( sink_flag && dist < R_SINK ){
            double t_visc = TAU_SINK;
            
//========================== Set t_visc to local viscous time ====================================
            if( sink_flag==2 ){
                double nu = explicit_viscosity;
                double M = pl->M;
                double eps = G_EPS;
                if( alpha_flag ){
                    double alpha = explicit_viscosity;
                    double c = sqrt( gamma_law*prim[PPP]/prim[RHO] );
                    double h = c*pow( r , 1.5 );
                    //viscosity calculation for Binary system--based on Yike's
                    if( cs2Choice==4 ){
                        //Ignoring contribution from other BH b/c in sink 
                        double denom = 1./sqrt( M*pow(dist*dist+eps*eps,-1.5) ); 
                        h = c*denom;
                    }
                    nu = alpha*c*h;
                }          
                t_visc = 2./3.*dist*dist/nu;  //TODO:Might need a factor of u^-0.5 here...
                //printf("r = %f, tau = %f\n", r, t_visc );
            }
//=================================================================================================
            if( t_visc < 10.*dt )
                t_visc = 10.*dt;
            double sink_frac = dt/t_visc;
        
            ////Track mass that will be eaten by sink
            pl->m_accr += sink_frac*cons[DDD];
            
            //printf("Mass removed, Mdot(t): (%e, %e)\t [%e, %e, %e]\n", sink_frac*cons[DDD], cons[DDD]/t_visc, cons[DDD]/dV, dV, sink_frac);

            
            //prim[RHO] -= sink_frac*rho;
            cons[DDD] -= sink_frac*cons[DDD];
            cons[SRR] -= sink_frac*cons[SRR];
            cons[LLL] -= sink_frac*cons[LLL]*pefficiency;
            cons[SZZ] -= sink_frac*cons[SZZ];
            cons[TAU] -= sink_frac*cons[TAU];
            //int q;
            //double floor = 1e-5;
            //for( q=0; q<NUM_Q; q++ ){i
            //  cons[q] -= sink_frac*cons[q];
            //}
            
            //Calculate Accretion Torques
            double mdot = cons[DDD]/t_visc;  //prim[RHO]/t_visc; //cons[DDD] = prim[RHO]*dV
            torques[1] += mdot*rbh*(r*prim[UPP] - vbh); //sign chose to match yike... (correct?)
        }
        double fr, fp, fz;
        planetaryForce( pl, r, phi, 0.0, &fr, &fp, &fz, 1);
        //fr, fp and fz are on fluid element, so force on bh is -fj
        torques[0] += cons[DDD]*(-fp)*rbh; //*dV? //(rho-1.0)???
        if( abs(cons[DDD]*fp*rbh) > 100 )
            printf("Big torque (%e) at r = %e", cons[DDD]*(-fp)*rbh, r);
    }
}

void density_sink_yike_old( struct domain *theDomain, double *prim, double *cons, double *xp, double *xm, double dV, double dt ){

    double rho = prim[RHO];
    double vr  = prim[URR];
    double om  = prim[UPP];
    double vz  = prim[UZZ];
    double rp = xp[0];
    double rm = xm[0];
    //double r = get_moment_arm( xp, xm );
    double r = 0.5*(xp[0]+xm[0]);
    double dphi = get_dp( xp[1],xm[1] );
    double phip = xp[1];
    double phim = phip - dphi;
    double phi  = 0.5*(phip+phim);
    double dist = nearest_planet_dist( theDomain, 0.5*(rp+rm), 0.5*(phip+phim) );

    if( sink_flag && dist < R_SINK ){
        double t_visc = TAU_SINK;
        if( t_visc < 10.*dt )
            t_visc = 10.*dt;
        double sink_frac = dt/t_visc;
        
        //prim[RHO] -= sink_frac*rho;
        cons[DDD] -= sink_frac*cons[DDD];
        cons[SRR] -= sink_frac*cons[SRR];
        cons[LLL] -= sink_frac*cons[LLL];
        cons[SZZ] -= sink_frac*cons[SZZ];
        cons[TAU] -= sink_frac*cons[TAU];
        //int q;
        //double floor = 1e-5;
        //for( q=0; q<NUM_Q; q++ ){i
        //  cons[q] -= sink_frac*cons[q];
        //}

    }
}

void initial( double*, double* );

void damping( double *cons, double *xp, double *xm, double dV, double dt ){
   
    //TODO: Maybe put everything in the if statement so don't do it if don't have to
    //TODO: Track how much ang mom we're damping against e.g. what is (x-x0)*damp_frac?
    
    //Make a prims array and call prim to cons
    double prim0[NUM_Q];
    double cons0[NUM_Q];

    double rp = xp[0];
    double rm = xm[0];
    double r = get_moment_arm(xp,xm);
    double dphi = get_dp( xp[1], xm[1] );
    double phip = xp[1];
    double phim = phip-dphi;

    //Set Prim0's to initial conditions
    double x[3] = { r, 0.5*(phip+phim), 0.5*(xp[2]+xm[2]) };
    initial( prim0, x );

    //Build Cons0
    prim2cons( prim0, cons0, x, dV );

    //set timescale tau
    double OM_P = 2*M_PI;
    double tau = 100.0*OM_P;
    //double tau = 2*M_PI*pow(r_max, 1.5);
    if( tau < 10*dt )
        tau = 10*dt;

    double damp_frac = dt/tau;
    if( r_damp>0.0 && r>r_damp ){
        //do the damping!
        int q;
        for( q=0; q<NUM_Q; ++q ){
            cons[q] -= damp_frac*(cons[q] - cons0[q])*Theta(r);
        }
    }
}



void visc_flux( double * prim , double * gprim , double * flux , double * x , double * n, struct planet *thePlanets ){

   double r = x[0];
   double nu = explicit_viscosity;

//================ Yike Binary Visc. Calc. =================
//TODO:   For eccentric or live binaries, need to fix this
//          -> source() has acces to planet objects
/*
   double q = q_planet;
   double om_bh = pow( a, -1.5 );
   double mu = q/(1.+q);
   double eps = G_EPS;

   double M1 = 1.-mu;
   double r1 = a*mu;
   double phi1 = om_bh*t_step - M_PI;
   double d1 = get_planet_2dist( r, x[1], r1, phi1 );

   double M2 = mu;
   double r2 = a*(1.-mu);
   double phi2 = om_bh*t_step;
   double d2 = get_planet_2dist( r, x[1], r2, phi2 );
*/
   double eps = G_EPS;

   double M1 = thePlanets[0].M;
   double M2 = thePlanets[1].M;
   double d1 = get_planet_2dist( r, x[1], thePlanets[0].r, thePlanets[0].phi );
   double d2 = get_planet_2dist( r, x[1], thePlanets[1].r, thePlanets[1].phi );
//==========================================================

   if( alpha_flag ){
        double alpha = explicit_viscosity;
        double c = sqrt( gamma_law*prim[PPP]/prim[RHO] );
        double h = c*pow( r , 1.5 );
        //viscosity calculation for Binary system--based on Yike's
        if( cs2Choice==4 ){
           double denom = 1./sqrt( M1*pow(d1*d1+eps*eps,-1.5) + M2*pow(d2*d2+eps*eps,-1.5) );
           h = c*denom;
        }
        nu = alpha*c*h;
    }
/*
   if( alpha_flag ){
      double alpha = explicit_viscosity;
      double c = sqrt( gamma_law*prim[PPP]/prim[RHO] );
      double h = c*pow( r , 1.5 );
      
      nu = alpha*c*h;
   }
*/

   double rho = prim[RHO];
   double vr  = prim[URR];
   double om  = prim[UPP];
   double om_off = om - get_om(r);
   double vz  = prim[UZZ];

   double dnvr = gprim[URR];
   double dnom = gprim[UPP];
   double dnvz = gprim[UZZ];

   flux[SRR] = -nu*rho*( dnvr - n[1]*2.*om );
   flux[LLL] = -nu*rho*( r*r*dnom + n[1]*2.*vr );
   flux[SZZ] = -nu*rho*dnvz;
   flux[TAU] = -nu*rho*( vr*dnvr+r*r*om_off*dnom+vz*dnvz );//- 2.*r*om_off*om );

}

void flux_to_E( double * Flux , double * Ustr , double * x , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann , int dim ){

   //Silence is Golden.

}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss , double * n , double * x , double * Bpack ){

   double r = x[0];
   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[URR]*n[0] + prim1[UPP]*n[1]*r + prim1[UZZ]*n[2];

   double cs1 = sqrt(gamma_law*P1/rho1);

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[URR]*n[0] + prim2[UPP]*n[1]*r + prim2[UZZ]*n[2];

   double cs2 = sqrt(gamma_law*P2/rho2);

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;

}

double get_dL( double * , double * , int );

double mindt(double * prim , double w , double * xp , double * xm ){

   double r = .5*(xp[0]+xm[0]);
   double Pp  = prim[PPP];
   double rho = prim[RHO];
   double vp  = (prim[UPP]-w)*r;
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double cs  = sqrt(gamma_law*Pp/rho);

   
   double maxvr = cs + fabs(vr);
   double maxvp = cs + fabs(vp);
   double maxvz = cs + fabs(vz);

   double dtr = get_dL(xp,xm,1)/maxvr;
   double dtp = get_dL(xp,xm,0)/maxvp;
   double dtz = get_dL(xp,xm,2)/maxvz;
   
   double dt = dtr;
   if( dt > dtp ) dt = dtp;
   if( dt > dtz ) dt = dtz;
/*
   double dL0 = get_dL(xp,xm,0);
   double dL1 = get_dL(xp,xm,1);
   double dL2 = get_dL(xp,xm,2);
   double dx = dL0;
   if( dx>dL1 ) dx = dL1;
   if( dx>dL2 ) dx = dL2;

   double nu = explicit_viscosity;

   if( alpha_flag ){
      double alpha = explicit_viscosity;
      double c = sqrt( gamma_law*prim[PPP]/prim[RHO] );
      double h = c*pow( r , 1.5 );
      nu = alpha*c*h;
   }

   double dt_visc = .03*dx*dx/nu;
   if( dt > dt_visc ) dt = dt_visc;
*/
   double delta_t = 1e-4;
   if( dt <= delta_t ){
        //printf("rho = %f\n", rho);
        //printf("P   = %f\n", Pp );
        //printf("cs  = %f\n", cs );
        //printf("Fluid velocity in cell: (%f, %f, %f)\n", vr, vp, vz);
   }
   
   return( dt );

}
/*
double getReynolds( double * prim , double w , double * x , double dx ){

   double r = x[0];
   double nu = explicit_viscosity;

   double vr = prim[URR];
   double omega = prim[UPP];
   double vp = omega*r-w;
   double vz = prim[UZZ];

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double cs = sqrt(gamma_law*Pp/rho);
   
   double v = sqrt(vr*vr + vp*vp + vz*vz);

   double Re = (v+cs)*dx/nu;
   
   return(Re);

}
*/

void reflect_prims(double * prim, double * x, int dim)
{
    //dim == 0: r, dim == 1: p, dim == 2: z
    if(dim == 0)
        prim[URR] = -prim[URR];
    else if(dim == 1)
        prim[UPP] = -prim[UPP];
    else if(dim == 2)
        prim[UZZ] = -prim[UZZ];
}
