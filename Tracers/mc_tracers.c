#include "../paul.h"

//static int Nmc = 0.0;  //I think unnecessary now--Ntr becomes Nmc for mc tracers

double get_dV( double *, double * );

void setTracerParams( struct domain * theDomain){

   //int initType = theDomain->theParList.tr_init_type;
   //int mc_flag = theDomain->theParList.tr_init_type;
   //int Ncells = theDomain->Ncells;
   int Nmc = theDomain->theParList.num_tracers;
   theDomain->Ntr = Nmc;

   //theDomain->tr_out = theDomain->theParList.tr_out_flag;
   MPI_Barrier( theDomain->theComm );
   printf("Rank %d: Ntr = %d \n", theDomain->rank, theDomain->Ntr);
}


void printTracerCoords( struct domain * );
void syncTracerIDs( struct domain* );
double get_moment_arm( double*, double* );
void addTracer( struct tracerList * );

void init_tracerList( struct domain *theDomain ){
   //For MC tracers: Initializes tracer lists for each cell
   struct cell **theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   int Nmc = theDomain->Ntr;
   double *r_jph = theDomain->r_jph;
   double *z_kph = theDomain->z_kph;

   int id = 0;
   int i,j,k,n;
   for( j=0; j<Nr; ++j){
      for( k=0; k<Nz; ++k ){
         int jk = j+Nr*k;
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double zm = z_kph[k-1];
         double zp = z_kph[k];
         for( i=0; i<Np[jk]; ++i ){
            struct cell *c = &(theCells[jk][i]);
            c->myList = (struct tracerList *) malloc( sizeof(struct tracerList) );
            c->myList->head = NULL;
            c->myList->size = 0;

            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {rp,phip,zp};
            double xm[3] = {rm,phim,zm};
            double r = get_moment_arm(xp,xm);

            for( n=0; n<Nmc; ++n ){ //create local tracer list
               addTracer( c->myList );
               struct tracer *tr = c->myList->head;
               tr->R   = r; c->coords[0] = r;
               tr->Phi = 0.5*(phim+phip); c->coords[1] = tr->Phi;
               tr->Z   = 0.5*(zm+zp); c->coords[2] = tr->Z;
               tr->ID  = id;
               tr->Type = 1;
               tr->rmFlag = 0;
               //tr->myCell = c;   //don't think I need this anymore
               id++;
            }
         }
      }
   }

//   syncTracerIDs( theDomain );  //number of tracers on the process is just Ncells*Nmc
}


void initializeTracers( struct domain *theDomain ){
/*
  int initType = theDomain->theParList.tr_init_type;
  int id = 0.0;
  int rank = theDomain->rank;
  struct tracer *tr = theDomain->theTracers->head;

   if( !rank )
      printf("Initializing tracers per cell\n");

   struct cell **theCells = theDomain->theCells;
   int  Nr = theDomain->Nr;
   int  Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   double *r_jph = theDomain->r_jph;
   double *z_kph = theDomain->z_kph;

   int i,j,k,n;
   for( j=0; j<Nr; ++j ){
     for( k=0; k<Nz; ++k ){
         int jk = j+Nr*k;
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double zm = z_kph[k-1];
         double zp = z_kph[k];
         for( i=0; i<Np[jk]; ++i ){
            struct cell *c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {rp,phip,zp};
            double xm[3] = {rm,phim,zm};
            double r = get_moment_arm(xp,xm);
            for( n=0; n<Nmc; ++n ){
               if( tr==NULL ){
                  printf("Tracer list size does not match (cell number)*(tracers per cell)\n");
                  break;
               }   
               tr->R   = r;
               tr->Phi = 0.5*(phim+phip);
               tr->Z   = 0.5*(zm+zp);
               tr->ID  = id;
               tr->Type   = 1;
               tr->rmFlag = 0;
               tr->myCell = c;   //Not sure I'm doing the pointer to the cell right..
               id++;
               tr = tr->next;
            }
         }
      }
   }
 //printTracerCoords( theDomain );
  syncTracerIDs( theDomain );
*/
}

void transferTracer( struct tracer *tr,  struct cell *c ){
   
   struct tracerList *list = c->myList;
   struct tracer *temp = list->head;
   list->head = tr;
   tr->R = c->coords[0];
   tr->Phi = c->coords[1];
   tr->Z = c->coords[2];
   tr->next = temp;
   list->size += 1;
   printf("Tracer was moved!!!\n");

}

void tracerMC( double mflux, struct cell *cL, struct cell *cR ){
   
   double M_i  = cL->prim[0];
   double delM = mflux;
   double p_flux = delM/M_i; //need to check if this is positive or negative

   double x;
   struct tracerList *list = cL->myList;
   struct tracer *tr = list->head, *prev=NULL, *temp=NULL;
   if( p_flux > 0  ){  //if outgoing mass flux
      while( tr!=NULL ){
         x = rand()/RAND_MAX;
         temp = tr->next;
         if( x<p_flux ){
            if( prev )
               prev->next = tr->next;
            else
               list->head = tr->next;
            transferTracer( tr, cR );
            list->size -= 1;
         }
         else{
            prev = tr;
         }
         tr = temp;
      }
   }

}


void updateTracers(struct domain *theDomain, double dt){
   return;
}

void tracer_RK_copy( struct tracer * tr){
   return;
}

void tracer_RK_adjust( struct tracer *tr, double RK ){
   return;
}


int trOutStep( struct domain *theDomain ){

   int step = theDomain->mdStep;
   int printStep = theDomain->theParList.tr_out_step;

   if( (step+printStep)%printStep==0 )
      return 1;
   else
      return 0;
}

int getListSize( struct tracerList * );

void tracerOutput( struct domain *theDomain ){

   char filename[256];
   sprintf(filename, "%s.xyz", "tracers_mc" );

   struct cell **theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;

   int Ntr_tot = theDomain->Ntr*theDomain->Ncells;
   int step = theDomain->mdStep;
   int rank = theDomain->rank;
   int size = theDomain->size;
   double t = theDomain->t;

   int outFlag = theDomain->theParList.tr_out_flag;

   MPI_Allreduce( MPI_IN_PLACE, &Ntr_tot, 1, MPI_INT, MPI_SUM, theDomain->theComm );

   if( rank==0 && step==0 ){
	   FILE * pFile = fopen(filename, "w");
      printf("Rewriting Tracer Out-File\n");
	   fclose(pFile);
   }
   MPI_Barrier( theDomain->theComm);

   int rk;
   int count=0;
   for( rk=0; rk<size; ++rk){
      //MPI_Barrier( theDomain->theComm );
      if( rank==rk ){
         FILE * pFile = fopen(filename, "a");
         if( rank==0 && outFlag==0 ){
            fprintf(pFile, "%d \nAtoms. Timestep: %d\n", Ntr_tot+1, step);
            //printf(" Ntr_tot from output: %d \n", Ntr_tot );
            // step,time, id,tpe,  x, y, z, r, phi, vr, om, vz
            fprintf(pFile, "%d %4.2f %d %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n",
                     step,t, 0,0, 0.0,0.0,0.0, 0.0,0.0 ,0.0,0.0,0.0);
         }
         int i,j,k,jk;
         for( j=0; j<Nr; ++j){
            for( k=0; k<Nz; ++k){
               jk = j + k*Nr;
               for( i=0; i<Np[jk]; ++i){
                  struct cell *c = &(theCells[jk][i]);
                  struct tracer *tr = c->myList->head;
                  while( tr!=NULL ){
                     int id   = tr->ID;
   	               int type = tr->Type;
   	               double r = tr->R;
               	   double phi = tr->Phi;
      	            double z = tr->Z;
               	   double x = r*cos(phi);
                  	double y = r*sin(phi);
                  	double vr = tr->Vr;
      	            double om = tr->Omega;
                  	double vz = tr->Vz;
      	            fprintf(pFile, "%d %4.2f %d %d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n",                            step,t, 
                           id,type, 
                           x,y,z, 
                           r,phi, 
                           vr,om,vz);
                     tr = tr->next;
                     count++;

                  }
               }
            }
         }
         fclose(pFile);
      }
      MPI_Barrier( theDomain->theComm );
   }

}




