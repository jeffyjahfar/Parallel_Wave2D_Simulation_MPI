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

void printCurrentState(double *zt, int size){
	for(int i=0;i<size;i++){
	  for(int j=0;j<size;j++){
		  cout<<zt[i*size + j]<<" ";
	  }
	  cout<<"\n";
	}
	cout<<"\n";
}

//method to initialize wave for t=0
void initializeForTimeZero(double *zt0,int stripe, int size, int weight, int rank){
//	#pragma omp parallel for default( none ) shared( zt0,stripe,size,weight,rank )
	#pragma omp parallel for
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

//method to initialize wave for t=1 for blocks inside the stripe excluding border
void initializeForTimeOneInternal(double *zt1,double *zt0,int stripe, int size, double w1, int rank){
//	#pragma omp parallel for default( none ) shared( zt1,zt0,stripe,size,w1,rank )
	#pragma omp parallel for
	for( int i = 1; i < stripe-1; i++ ) {
	  for( int j = 0; j < size; j++ ) {
		int i_pos = rank*stripe+i;
		if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
			zt1[i*size + j] = 0.0;
		} else {
			zt1[i*size + j] = zt0[i*size + j]+ w1*(zt0[(i+1)*size + j]+zt0[(i-1)*size + j]+zt0[i*size + j+1]+zt0[i*size + j-1]-4*zt0[i*size + j]);
		}
	  }
	}
}

//method to initialize wave for t=1 for blocks in the border
void initializeForTimeOneBorders(double *zt1,double *zt0,double *ztl,double *ztr, int stripe, int size, double w1, int rank){
//	#pragma omp parallel for default( none ) shared( zt1,zt0,ztl,ztr,stripe,size,w1,rank )
	#pragma omp parallel for
	for( int i = 0; i < stripe; i+=stripe-1 ) {
	  for( int j = 0; j < size; j++ ) {
		int i_pos = rank*stripe+i;
		if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
			zt1[i*size + j] = 0.0;
		} else {
			if(i != 0){
				  zt1[i*size + j] = zt0[i*size + j]+ w1*(ztr[j]+zt0[(i-1)*size + j]+zt0[i*size + j+1]+zt0[i*size + j-1]-4*zt0[i*size + j]);
			}else{
				  zt1[i*size + j] = zt0[i*size + j]+ w1*(zt0[(i+1)*size + j]+ztl[j]+zt0[i*size + j+1]+zt0[i*size + j-1]-4*zt0[i*size + j]);
			}
		}
	  }
	}
}

//method to initialize wave for t=2 to max_time for blocks inside the stripe excluding border
void initializeForTimeTInternal(double *zt2, double *zt1, double *zt0, int stripe, int size, double w2, int rank){
//	#pragma omp parallel for default( none ) shared( zt2,zt1,zt0,stripe,size,w2,rank )
	#pragma omp parallel for
	for( int i = 1; i < stripe-1; i++ ) {
		for( int j = 0; j < size; j++ ) {
			int i_pos = rank*stripe+i;
			if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
				zt2[i*size + j] = 0.0;
			} else {
				zt2[i*size + j] = 2*zt1[i*size + j]-zt0[i*size + j]+w2*(zt1[(i+1)*size + j]+zt1[(i-1)*size + j]+zt1[i*size + j+1]+zt1[i*size + j-1]-4*zt1[i*size + j]);
			}
		}
	}
}

//method to initialize wave for t=2 to max_time for blocks in the border
void initializeForTimeTBorders(double *zt2, double *zt1,double *ztl, double *ztr, double *zt0, int stripe, int size, double w2, int rank){
//	#pragma omp parallel for default( none ) shared( zt2,zt1,ztl,ztr,zt0,stripe,size,w2,rank )
	#pragma omp parallel for
	for( int i = 0; i < stripe; i+=stripe-1 ) {
		for( int j = 0; j < size; j++ ) {
			int i_pos = rank*stripe+i;
			if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
				zt2[i*size + j] = 0.0;
			} else {
				if(i != 0){
					zt2[i*size + j] = 2*zt1[i*size + j]-zt0[i*size + j]+w2*(ztr[j]+zt1[(i-1)*size + j]+zt1[i*size + j+1]+zt1[i*size + j-1]-4*zt1[i*size + j]);
				}else{
					zt2[i*size + j] = 2*zt1[i*size + j]-zt0[i*size + j]+w2*(zt1[(i+1)*size + j]+ztl[j]+zt1[i*size + j+1]+zt1[i*size + j-1]-4*zt1[i*size + j]);
				}
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

	// change # of threads
	omp_set_num_threads( nThreads );

	#pragma omp parallel for default( none ) shared( z,size )
	for ( int p = 0; p < 3; p++ ) {
	  for ( int i = 0; i < size; i++ ){
		  for ( int j = 0; j < size; j++ ){
			  z[p][i][j] = 0.0; // no wave
		  }
	  }
	}

	// start MPI
	MPI_Init( &argc, &argv );
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


  // time = 0;
  // initialize the simulation space: calculate z[0][][]
  int weight = size / default_size;
  MPI_Bcast( &size, 1, MPI_INT, 0, MPI_COMM_WORLD );

  int stripe = size/mpi_size;

	if ( my_rank == 0 ) {
		cerr << "rank["<<my_rank<<"]'s range = "<<stripe*my_rank<<" ~ "<<stripe*(my_rank+1)-1<< endl;
//		#pragma omp parallel for
		for (int rank=1; rank<mpi_size; rank++){
			cerr << "rank["<<rank<<"]'s range = "<<stripe*rank<<" ~ "<<stripe*(rank+1)-1<< endl;
			MPI_Send(*zt0+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,3, MPI_COMM_WORLD);
			MPI_Send(*zt1+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,4, MPI_COMM_WORLD);
			MPI_Send(*zt2+size*stripe*rank, stripe*size, MPI_DOUBLE,rank,5, MPI_COMM_WORLD);
		}
	}else{
		int source = 0;
		MPI_Recv(zt0, stripe*size, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);
		MPI_Recv(zt1, stripe*size, MPI_DOUBLE, source, 4, MPI_COMM_WORLD, &status);
		MPI_Recv(zt1, stripe*size, MPI_DOUBLE, source, 5, MPI_COMM_WORLD, &status);
	}

  initializeForTimeZero(*zt0, stripe, size, weight,my_rank);

//  if(my_rank == 0){
////	  #pragma omp parallel for
//	  for (int rank=1; rank<mpi_size; rank++){
//		MPI_Recv(*zt0+size*stripe*rank, stripe*size, MPI_DOUBLE, rank, 2, MPI_COMM_WORLD, &status);
//	  }
////	  cout<<0<<"\n";
////	  printCurrentState(*zt0, size);
//  }else{
//	  MPI_Send(zt0, stripe*size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
//  }

  // time = 1
  // calculate z[1][][] 

  //Initialize constants in the equation
  double w1 = (c*c/2)*(dt*dt)/(dd*dd);

  //Loop to initialize all blocks for t=1
  initializeForTimeOneInternal(*zt1, *zt0, stripe, size, w1,my_rank);

	double *ztl = new double[1*size];
	double *ztr = new double[1*size];
  	if(my_rank %2==0){
  		//send left
  			if(my_rank !=0){
  				MPI_Send(*zt0+1*size, 1*size, MPI_DOUBLE, my_rank-1, 6, MPI_COMM_WORLD);
  			}

  			//send right
  			if(my_rank != mpi_size-1){
  				MPI_Send(*zt0+(stripe-1)*size, 1*size, MPI_DOUBLE, my_rank+1, 7, MPI_COMM_WORLD);
  			}
  			//receive left
  			if(my_rank !=0){
				MPI_Recv(ztl, 1*size, MPI_DOUBLE, my_rank-1, 7, MPI_COMM_WORLD, &status);
  			}

  			//receive right
  			if(my_rank !=mpi_size-1){
  				MPI_Recv(ztr, 1*size, MPI_DOUBLE, my_rank+1, 6, MPI_COMM_WORLD, &status);
  			}
  	}else{

			//receive left
			if(my_rank !=0){
				MPI_Recv(ztl, 1*size, MPI_DOUBLE, my_rank-1, 7, MPI_COMM_WORLD, &status);
			}

  			//receive right
  			if(my_rank !=mpi_size-1){
  				MPI_Recv(ztr, 1*size, MPI_DOUBLE, my_rank+1, 6, MPI_COMM_WORLD, &status);
  			}

  			//send left
			if(my_rank !=0){
				MPI_Send(*zt0+1*size, 1*size, MPI_DOUBLE, my_rank-1, 6, MPI_COMM_WORLD);
			}

			//send right
			if(my_rank != mpi_size-1){
				MPI_Send(*zt0+(stripe-1)*size, 1*size, MPI_DOUBLE, my_rank+1, 7, MPI_COMM_WORLD);
			}
  	}

	initializeForTimeOneBorders(*zt1, *zt0,ztl, ztr, stripe, size, w1,my_rank);

//  if(my_rank == 0){
////	  #pragma omp parallel for
//	  for (int rank=1; rank<mpi_size; rank++){
//		int source = rank;
//		MPI_Recv(*zt1+size*stripe*rank, stripe*size, MPI_DOUBLE, source, 4, MPI_COMM_WORLD, &status);
//	  }
////	  cout<<1<<"\n";
////	  printCurrentState(*zt1, size);
//  }else{
//	  MPI_Send(zt1, stripe*size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
//  }

  // simulate wave diffusion from time = 2

  //Initialize constants in the equation
  double w2 = w1*2;

  //Loop to initialize all blocks for t=2 ....max_time-1
  for ( int t = 2; t <max_time; t++ ){

	  //initialize waves within the stripe
	  initializeForTimeTInternal(*zt2, *zt1, *zt0, stripe, size,w2,my_rank);

	  	  	if ( my_rank % 2 == 0 ) {
	  		  //send left
	  		  	if(my_rank !=0){
	  		  		MPI_Send(*zt1+1*size, 1*size, MPI_DOUBLE, my_rank-1, 6, MPI_COMM_WORLD);
	  		  	}

	  		  	//send right
	  		  	if(my_rank != mpi_size-1){
	  		  		MPI_Send(*zt1+(stripe-1)*size, 1*size, MPI_DOUBLE, my_rank+1, 7, MPI_COMM_WORLD);
	  		  	}

	  		  	//receive left
	  		  	if(my_rank !=0){
	  		  		MPI_Recv(ztl, 1*size, MPI_DOUBLE, my_rank-1, 8, MPI_COMM_WORLD, &status);
	  		  	}

	  		  	//receive right
	  		  	if(my_rank !=mpi_size-1){
	  		  		MPI_Recv(ztr, 1*size, MPI_DOUBLE, my_rank+1, 9, MPI_COMM_WORLD, &status);
	  		  	}
	  	  	}else {
	  		  	//receive left
	  		  	if(my_rank !=0){
	  		  		MPI_Recv(ztl, 1*size, MPI_DOUBLE, my_rank-1, 7, MPI_COMM_WORLD, &status);
	  		  	}

	  		  	//receive right
	  		  	if(my_rank !=mpi_size-1){
	  		  		MPI_Recv(ztr, 1*size, MPI_DOUBLE, my_rank+1, 6, MPI_COMM_WORLD, &status);
	  		  	}

	  		  //send left
				if(my_rank !=0){
					MPI_Send(*zt1+1*size, 1*size, MPI_DOUBLE, my_rank-1, 9, MPI_COMM_WORLD);
				}

				//send right
				if(my_rank != mpi_size-1){
					MPI_Send(*zt1+(stripe-1)*size, 1*size, MPI_DOUBLE, my_rank+1, 8, MPI_COMM_WORLD);
				}
	  	  	}

	  initializeForTimeTBorders(*zt2, *zt1,ztl,ztr, *zt0, stripe, size,w2,my_rank);

	  //aggregate all result at rank 0 for printing
	  if(interval != 0 && (t%interval==0 || t==max_time-1)){
		  if(my_rank == 0){
			  #pragma omp parallel for
			  for (int rank=1; rank<mpi_size; rank++){
				int source = rank;
				MPI_Recv(*zt2+size*stripe*rank, stripe*size, MPI_DOUBLE, source, 8, MPI_COMM_WORLD, &status);
			  }
			  //Print current state after gathering values from all ranks
			  cout<<t<<"\n";
			  printCurrentState(*zt2, size);

			}else{
			  MPI_Send(zt2, stripe*size, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
			}
	  }
	  //Rotate Pointers
	  double (*temp)[size] = zt0;
	  zt0 = zt1;
	  zt1 = zt2;
	  zt2 = temp;
  }
  // end of simulation

  // stop the timer
  if ( my_rank == 0 )
	  cerr << "Elapsed time = " << time.lap( ) << endl;

  MPI_Finalize( ); // shut down MPI

  return 0;
}
