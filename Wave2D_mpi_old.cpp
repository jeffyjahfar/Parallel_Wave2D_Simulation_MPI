#include <iostream>
#include <string.h>
#include "Timer.h"
#include <stdlib.h>   // atoi
#include "mpi.h"
#include <omp.h>     // OpenMP

int default_size = 100;  // the default system size
int defaultCellWidth = 8;
double c = 1.0;      // wave speed
double dt = 0.1;     // time quantum
double dd = 2.0;     // change in system

using namespace std;

//task parallelization
//void initializeWithZeroes(double *z,int stripe,int size){
//  for ( int p = 0; p < 3; p++ ) {
//	  for ( int i = 0; i < stripe; i++ ){
//		  for ( int j = 0; j < size; j++ ){
//			  (*z)[p][i][j] = 0.0; // no wave
//		  }
////		  memset(z[p][i],0,size);
//	  }
//  }
//}

void printCurrentState(double *zt, int size){
	for(int i=0;i<size;i++){
	  for(int j=0;j<size;j++){
		  cout<<zt[i*size + j]<<" ";
	  }
	  cout<<"\n";
	}
	cout<<"\n";
}

//task parallelization
void initializeForTimeZero(double *zt0,int stripe, int size, int weight, int rank){
	#pragma omp parallel for default( none ) shared( zt0,stripe,size,weight,rank )
	for( int i = 0; i < stripe; i++ ) {
		for( int j = 0; j < size; j++ ) {
			int i_pos = rank*stripe+i;
		  if( i_pos > 40 * weight && i_pos < 60 * weight  && j > 40 * weight && j < 60 * weight ) {
			  zt0[i * size + j] = 20.0;
		  } else {
			  zt0[i * size + j] = 0.0;
		  }
		}
	}
}

//task parallelization
void initializeForTimeOne(double *zt1,double *zt0,int stripe, int size, double w1, int rank){
	#pragma omp parallel for default( none ) shared( zt1,zt0,stripe,size,w1,rank )
	for( int i = 0; i < stripe; i++ ) {
	  for( int j = 0; j < size; j++ ) {
		int i_pos = i+rank*stripe;
		if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
			zt1[i*size + j] = 0.0;
		} else {
			zt1[i*size + j] = zt0[i_pos*size + j]+ w1*(zt0[(i_pos+1)*size + j]+zt0[(i_pos-1)*size + j]+zt0[i_pos*size + j+1]+zt0[i_pos*size + j-1]-4*zt0[i_pos*size + j]);
		}
	  }
	}
}

void initializeForTimeT(double *zt2, double *zt1, double *zt0, int stripe, int size, double w2, int rank){
	#pragma omp parallel for default( none ) shared( zt2,zt1,zt0,stripe,size,w2,rank )
	for( int i = 0; i < stripe; i++ ) {
		for( int j = 0; j < size; j++ ) {
			int i_pos = i+rank*stripe;
			if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
				zt2[i*size + j] = 0.0;
			} else {
				zt2[i*size + j] = 2*zt1[i_pos*size + j]-zt0[i*size + j]+w2*(zt1[(i_pos+1)*size + j]+zt1[(i_pos-1)*size + j]+zt1[i_pos*size + j+1]+zt1[i_pos*size + j-1]-4*zt1[i_pos*size + j]);
			}
		}
	}
}

int main( int argc, char *argv[] ) {
	int my_rank = 0;            // used by MPI
	int mpi_size = 1;           // used by MPI
	MPI_Status status;
	int tag = 0;              // tag for messages
	int nThreads = 1;

  // verify arguments
  if ( argc != 5 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    return -1;
  }
  int size = atoi( argv[1] );
  int max_time = atoi( argv[2] );
  int interval  = atoi( argv[3] );
  nThreads = atoi( argv[4] );

  if ( size < 100 || max_time < 3 || interval < 0 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
    return -1;
  }

	// create a simulation space
	double z[3][size][size];
	#pragma omp parallel for default( none ) shared( z,size )
	for ( int p = 0; p < 3; p++ ) {
	  for ( int i = 0; i < size; i++ ){
		  for ( int j = 0; j < size; j++ ){
			  z[p][i][j] = 0.0; // no wave
		  }
	//		  memset(z[p][i],0,size);
	  }
	}

	MPI_Init( &argc, &argv ); // start MPI
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );

	//Initialize pointers
	double (*zt0)[size] =  z[0];
	double (*zt1)[size] =  z[1];
	double (*zt2)[size] =  z[2];

  // start a timer
  Timer time;
  if(my_rank ==0)
	  time.start( );

  // change # of threads
    omp_set_num_threads( nThreads );

  // time = 0;
  // initialize the simulation space: calculate z[0][][]
  int weight = size / default_size;
  MPI_Bcast( &size, 1, MPI_INT, 0, MPI_COMM_WORLD );

  int stripe = size/mpi_size;

  if ( my_rank == 0 ) {
	  for (int rank=1; rank<mpi_size; rank++){
	        MPI_Send(*zt0+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,1, MPI_COMM_WORLD);
	  }
  }else{
	  MPI_Recv(zt0, stripe*size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
  }

  initializeForTimeZero(*zt0, stripe, size, weight,my_rank);

  if(my_rank == 0){
	  cerr << "rank["<<my_rank<<"]'s range = "<<stripe*my_rank<<" ~ "<<stripe*(my_rank+1)-1<< endl;
	  for (int rank=1; rank<mpi_size; rank++){
		cerr << "rank["<<rank<<"]'s range = "<<stripe*rank<<" ~ "<<stripe*(rank+1)-1<< endl;
		MPI_Recv(*zt0+size*stripe*rank, stripe*size, MPI_DOUBLE, rank, 2, MPI_COMM_WORLD, &status);
	  }
//	  cout<<0<<"\n";
//	  printCurrentState(*zt0, size);
  }else{
	  MPI_Send(zt0, stripe*size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
  }

  // time = 1
  // calculate z[1][][] 
  // cells not on edge
  // IMPLEMENT BY YOURSELF !!!
  //Initialize constants in the equation
  double w1 = (c*c/2)*(dt*dt)/(dd*dd);

  if ( my_rank == 0 ) {
	  for (int rank=1; rank<mpi_size; rank++){
	        MPI_Send(*zt1+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,3, MPI_COMM_WORLD);
	        MPI_Send(*zt0, size*size, MPI_DOUBLE, rank, 3, MPI_COMM_WORLD);
	  }

  }else{
	  int source = 0;
	  MPI_Recv(zt1, stripe*size, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);
	  MPI_Recv(zt0, size*size, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);
  }

  //Loop to initialize all blocks for t=1
  initializeForTimeOne(*zt1, *zt0, stripe, size, w1,my_rank);

  if(my_rank == 0){
	  for (int rank=1; rank<mpi_size; rank++){
		int source = rank;
		MPI_Recv(*zt1+size*stripe*rank, stripe*size, MPI_DOUBLE, source, 4, MPI_COMM_WORLD, &status);
	  }
//	  cout<<1<<"\n";
//	  printCurrentState(*zt1, size);
  }else{
	  MPI_Send(zt1, stripe*size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
  }


  // simulate wave diffusion from time = 2
  // IMPLEMENT BY YOURSELF !!!

  //Initialize constants in the equation
  double w2 = w1*2;

  //Loop to initialize all blocks for t=2 ....max_time-1
  for ( int t = 2; t <max_time; t++ ){

	  if ( my_rank == 0 ) {
		  for (int rank=1; rank<mpi_size; rank++){
			    MPI_Send(*zt0+size*stripe*rank, size*stripe, MPI_DOUBLE,rank,5, MPI_COMM_WORLD);
			    MPI_Send(*zt1, size*size, MPI_DOUBLE, rank, 6, MPI_COMM_WORLD);
		        MPI_Send(*zt2+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,7, MPI_COMM_WORLD);
		  }
	  }else{
		  int source = 0;
		  MPI_Recv(zt0, size*stripe, MPI_DOUBLE, source, 5, MPI_COMM_WORLD, &status);
		  MPI_Recv(zt1, size*size, MPI_DOUBLE, source, 6, MPI_COMM_WORLD, &status);
		  MPI_Recv(zt2, stripe*size, MPI_DOUBLE, source, 7, MPI_COMM_WORLD, &status);
	  }

	  initializeForTimeT(*zt2, *zt1, *zt0, stripe, size,w2,my_rank);

	  if(my_rank == 0){
		  for (int rank=1; rank<mpi_size; rank++){
			int source = rank;
			MPI_Recv(*zt2+size*stripe*rank, stripe*size, MPI_DOUBLE, source, 8, MPI_COMM_WORLD, &status);
		  }

		  //Print current state after gathering values from all ranks
		  if(interval != 0 && (t%interval==0 || t==max_time-1)){
	  		  cout<<t<<"\n";
	  		  printCurrentState(*zt2, size);
		  }

	  }else{
		  MPI_Send(zt2, stripe*size, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
	  }

	  //Rotate Pointers
	  zt0 = zt1;
	  zt1 = zt2;
	  zt2 = zt0;
  }
  // end of simulation

  // stop the timer
  if ( my_rank == 0 )
	  cerr << "Elapsed time = " << time.lap( ) << endl;

  MPI_Finalize( ); // shut down MPI

  return 0;
}
