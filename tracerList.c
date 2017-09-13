#include "paul.h"


void addTracer( struct tracerList *theList ){
   //adds tracer to the front of the list
   struct tracer *tr = (struct tracer *) malloc( sizeof(struct tracer) );

   tr->next = theList->head;

   theList->head = tr;
   theList->size += 1;
}

int getListSize( struct tracerList *theList ){

   int size = 0;
   struct tracer *tr = theList->head;
   while( tr!=NULL ){
      size++;
      tr = tr->next;
   }
  return size;
}

void rmTracers( struct domain *theDomain ){

   struct tracerList *theList = theDomain->theTracers;
   struct tracer *tr = theList->head;
   struct tracer *prev = NULL;
   struct tracer *del  = NULL;
   while( tr!=NULL ){
      if( tr->rmFlag ){
         if( prev ){ //if prev has a valid value
            prev->next = tr->next;
         }
         else{ //if prev is NULL, need to repoint head
            theList->head = tr->next;
         }
         del = tr;
         tr = tr->next;
         int type = tr->Type;
         free(del);
      }
      else{
         prev = tr;
         tr = tr->next;
      }
   }
}

void init_tracerList( struct domain *theDomain ){

   int Ntr = theDomain->Ntr;
   int Nmc = theDomain->Nmc;
   struct tracerList *theList = theDomain->theTracers;
   theList->head = NULL;
   theList->size = 0;

   int N = Ntr;
   if( Nmc ){
      printf("Using %g Monte Carlo Tracers per Cell", Nmc)
      N = Ntr*Nmc;
   }

   int n;
   for( n=0; n<N; ++n){     //create linked list of Ntr tracers
      addTracer( theList );
   }

   //printf("Trying to get list size (%d were added)\n", theList->size);
   //int size = getListSize( theList );
   //printf("List Size: %d\n", size);

}

void printTracerCoords( struct domain *theDomain ){

   struct tracer *tr = theDomain->theTracers->head;
   int type;
   double r, z, phi;
   printf("Type, (r, phi, z)\n");
   while ( tr!=NULL ){
      type = tr->Type;
      r = tr->R; z = tr->Z; phi = tr->Phi;
      printf("%d,  (%4.2f, %4.2f, %4.2f)\n", type, r, phi, z);
      tr = tr->next;
   }
}

void printTracerVels( struct domain *theDomain ){

   struct tracer *tr = theDomain->theTracers->head;
   int type;
   double vr, vz, vp;
   printf("Vels: (vr, om, vz)\n");
   while( tr!=NULL ){
      vr = tr->Vr; vp = tr->Omega; vz = tr->Vz;
      printf("(%4.2f, %4.2f, %4.2f)\n", vr,vp,vz);
      tr = tr->next;
   }
}
