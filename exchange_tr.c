#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "paul.h"

void printSendMsg( int srank, int rrank, struct tracer *tr ){

   printf("Rank %d trying to send tracer to rank %d. ", srank,rrank);
   printf("X: ( %4.2f, %4.2f, %4.2f )  V: ( %4.2f, %4.2f, %4.2f )\n",
            tr->R,tr->Phi,tr->Z, tr->Vr,tr->Omega,tr->Vz);
}

void printRecvMsg( int rank, struct tracer *tr ){
   printf("Rank %d received tracer:  ", rank);
   printf("X: ( %4.2f, %4.2f, %4.2f )  V: ( %4.2f, %4.2f, %4.2f )\n",
            tr->R,tr->Phi,tr->Z, tr->Vr,tr->Omega,tr->Vz);
}

int testTracerCoords( struct domain *theDomain ){

   int rank = theDomain->rank;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      double r   = tr->R;
      double phi = tr->Phi;
      double z   = tr->Z;
      printf("%d, %4.2f, %4.2f, %4.2f\n", rank,r,phi,z);

      if( isnan(r) && isnan(z) && isnan(phi) ){
         printf("Tracer was corrupted on rank %d\n", rank);
         //if(tr->Type==50) printf(" - It was passed\n");
         return 1;
      }
      tr = tr->next;
   }
   return 0;
}

struct tr_lite{

   int id;
   int type;
   double x[3];
   double v[3];

};

void generate_mpi_tr( MPI_Datatype *tr_mpi ){

   struct tr_lite test;
   int count = 4;
   int blocksize[]      = { 1, 1, 3, 3 };
   MPI_Datatype types[] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
   MPI_Aint offsets[4];

   offsets[0] = (char *)&(test.id) - (char *)(&test);
   offsets[1] = (char *)&(test.type) - (char *)(&test);
   offsets[2] = (char *)&(test.x)    - (char *)(&test);
   offsets[3] = (char *)&(test.v)    - (char *)(&test);

   MPI_Type_create_struct( count, blocksize, offsets, types, tr_mpi );
   MPI_Type_commit( tr_mpi );
}

void copy_tr_to_lite( struct tracer *tr, struct tr_lite *trl ){

   double x[3] = { tr->R , tr->Phi  , tr->Z  };
   double v[3] = { tr->Vr, tr->Omega, tr->Vz };
   memcpy( trl->x, x, 3*sizeof(double) );
   memcpy( trl->v, v, 3*sizeof(double) );
   trl->type = tr->Type;
   trl->id   = tr->ID;

}

void copy_lite_to_tr( struct tr_lite *trl, struct tracer *tr ){

   tr->ID    = trl->id;
   tr->Type  = trl->type;
   tr->R     = trl->x[0];
   tr->Vr    = trl->v[0];
   tr->Phi   = trl->x[1];
   tr->Omega = trl->v[1];
   tr->Z     = trl->x[2];
   tr->Vz    = trl->v[2];

}

void generate_sendbuffer_tr( struct domain *theDomain, int numL, int numR,  struct tr_lite *trL, struct tr_lite *trR, int dim ){

   int iL = 0;
   int iR = 0;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      if( tr->rmFlag==1 ){
         copy_tr_to_lite( tr, trL+iL );
         iL++;
      }
      if( tr->rmFlag==2 ){
         copy_tr_to_lite( tr, trR+iR );
         iR++;
      }
      tr = tr->next;
   }

}

void generate_tr_number( struct domain *theDomain, int *numL, int *numR, int dim ){

   double r0   = theDomain->r0;
   double z0   = theDomain->z0;
   double delr = theDomain->delr;
   double delz = theDomain->delz;
   double x0[2]   = { r0  , z0  };
   double delx[2] = { delr, delz};

   int iL = 0;
   int iR = 0;
   double coord = 1e5;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      if( dim==0 )
         coord = tr->R;
      if( dim==1 )
         coord = tr->Z;
      if( coord < x0[dim] ){
         tr->rmFlag = 1;
         iL++;
         //printSendMsg( theDomain->rank, theDomain->left_rank[dim], tr );
         //printf("Sending tracer in dimension %d\n", dim);
      }
      if( x0[dim]+delx[dim] < coord ){
         tr->rmFlag = 2;
         iR++;
         //printSendMsg( theDomain->rank, theDomain->right_rank[dim], tr );
         //printf("Sending tracer in dimension %d\n", dim);
      }
      tr = tr->next;
   }
   *numL = iL;
   *numR = iR;

}

void addTracer( struct tracerList * );
void tracer_RK_copy( struct tracer * );
int getListSize( struct tracerList * );

void fill_tr_list( struct tracerList *theList, int ntr, struct tr_lite *tr_recv, int rank ){

   int n;
   struct tracer *tr = theList->head;
   for( n=0; n<ntr ; ++n){
      addTracer( theList );
      tr = theList->head;
      copy_lite_to_tr( (tr_recv + n) , tr );
      tracer_RK_copy( tr );
      tr->rmFlag = 0;
      tr->myCell = NULL; //do a new cell search and find its new cell?....
      //printRecvMsg( rank, tr );
   }

}

void rmTracers( struct domain * );
void printTracerCoords( struct domain * );

void exchangeTracers( struct domain *theDomain, int dim ){

   MPI_Datatype tr_mpi = {0};
   generate_mpi_tr( &tr_mpi );

   MPI_Comm grid_comm = theDomain->theComm;
   int  rank = theDomain->rank;
   int *left_rank = theDomain->left_rank;
   int *right_rank = theDomain->right_rank;

   int tag = 0;
   MPI_Status status;

   int numL, numR;
   int send_sizeL = 0;
   int send_sizeR = 0;
   int recv_sizeL = 0;
   int recv_sizeR = 0;

//------Count Tracers to Send----------------------------
   generate_tr_number( theDomain, &numL, &numR, dim );
   send_sizeL = numL;
   send_sizeR = numR;

//------Tell Neighbors how many to expect----------------
   MPI_Sendrecv( &send_sizeL, 1, MPI_INT,  left_rank[dim], tag,
                 &recv_sizeR, 1, MPI_INT, right_rank[dim], tag,
                 grid_comm, &status );
   MPI_Sendrecv( &send_sizeR, 1, MPI_INT, right_rank[dim], tag+1,
                 &recv_sizeL, 1, MPI_INT,  left_rank[dim], tag+1,
                 grid_comm, &status );

//------Build list of tracers to send---------------------
   struct tr_lite trL_send[send_sizeL];
   struct tr_lite trR_send[send_sizeR];
   struct tr_lite trL_recv[recv_sizeL];
   struct tr_lite trR_recv[recv_sizeR];
   generate_sendbuffer_tr( theDomain, numL, numR, trL_send, trR_send, dim );

//------Send lists of tracers-----------------------------
   MPI_Sendrecv( trL_send, send_sizeL, tr_mpi,  left_rank[dim], tag+2,
                 trR_recv, recv_sizeR, tr_mpi, right_rank[dim], tag+2,
                 grid_comm, &status );
   MPI_Sendrecv( trR_send, send_sizeR, tr_mpi, right_rank[dim], tag+3,
                 trL_recv, recv_sizeL, tr_mpi,  left_rank[dim], tag+3,
                 grid_comm, &status );
//------Add new tracers recieved--------------------------

/*   if( recv_sizeL!=0 || recv_sizeR !=0 ){
      printf("Before particles added\n Rank: %d\n", rank);
      printTracerCoords( theDomain );
   }
*/

   fill_tr_list( theDomain->theTracers, recv_sizeR, trR_recv, rank );
   fill_tr_list( theDomain->theTracers, recv_sizeL, trL_recv, rank );
   
/*   
   if( recv_sizeL!=0 || recv_sizeR!=0 ){
      printf("After tracers added\n");
      printTracerCoords( theDomain );
   }
*/

//------Remove tracers sent away--------------------------
   rmTracers( theDomain );

/*
   if( recv_sizeL!=0 || recv_sizeR!=0 ){
      printf("After tracers Removed\n");
      printTracerCoords( theDomain );
   }
*/
  
   MPI_Type_free( &tr_mpi );

/*
   if( recv_sizeL!=0 || recv_sizeR!=0 ){
      printf("After TypeFree\n");
      printTracerCoords( theDomain );
   }
   int rk;
   int brake = 0;
   for( rk=0; rk<size; rk++ ){
      MPI_Barrier(grid_comm);
      if( rk==rank )
         brake = testTracerCoords( theDomain );
   }
   if( brake ){
      //printf("Corruption happend in removing from the list\n");
      exit(0);
   }
*/
}

void syncTracerIDs( struct domain *theDomain ){

   MPI_Comm theComm = theDomain->theComm;
   int rank = theDomain->rank;
   int size = theDomain->size;
   int Ntr  = theDomain->Ntr;
   int tr_nums[size];

   //Let each rank know how many tracers are on all other ranks
   tr_nums[rank] = Ntr;
   MPI_Allgather( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 
                  tr_nums, 1, MPI_INT,
                  theComm );

   //Calculate number of tracers before your rank
   int i, sum=0;
   for( i=0; i<rank; i++ )
      sum += tr_nums[i];
   
   //Fix tracer IDs
   int id;
   struct tracer *tr = theDomain->theTracers->head;
   while( tr!=NULL ){
      id =  tr->ID;
      id += sum;
      tr->ID = id;
      tr = tr->next;
   }

/*   
   printf("Rank %d, Ntr: %d \n", rank,Ntr );
   if(rank==0){
      for( i=0; i<size; i++)
         printf("%d ", tr_nums[i]);
      printf("\n");
   }
*/
   //free(&tr_nums);

}
