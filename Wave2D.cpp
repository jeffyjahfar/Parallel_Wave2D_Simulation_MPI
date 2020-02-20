#include <iostream>
#include "Timer.h"
#include <stdlib.h>   // atoi

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

//task parallelization
void initializeForTimeZero(double *zt0, int size, int weight){
  for( int i = 0; i < size; i++ ) {
	for( int j = 0; j < size; j++ ) {
		int i_pos = i;
	  if( i_pos > 40 * weight && i_pos < 60 * weight  && j > 40 * weight && j < 60 * weight ) {
		  zt0[i * size + j] = 20.0;
	  } else {
		  zt0[i * size + j] = 0.0;
	  }
	}
  }
}

//task parallelization
void initializeForTimeOne(double *zt1,double *zt0, int size, double w1){
	for( int i = 0; i < size; i++ ) {
	  for( int j = 0; j < size; j++ ) {
		int i_pos = i;
		if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
			zt1[i*size + j] = 0.0;
		} else {
			zt1[i*size + j] = zt0[i_pos*size + j]+ w1*(zt0[(i_pos+1)*size + j]+zt0[(i_pos-1)*size + j]+zt0[i_pos*size + j+1]+zt0[i_pos*size + j-1]-4*zt0[i_pos*size + j]);
		}
	  }
	}
}

void initializeForTimeT(double *zt2, double *zt1, double *zt0, int size, double w2){
	for( int i = 0; i < size; i++ ) {
		for( int j = 0; j < size; j++ ) {
			int i_pos = i;
			if( i_pos==0 || j==0 || i_pos==size-1 || j==size-1 ) {
				zt2[i*size + j] = 0.0;
			} else {
				zt2[i*size + j] = 2*zt1[i_pos*size + j]-zt0[i_pos*size + j]+w2*(zt1[(i_pos+1)*size + j]+zt1[(i_pos-1)*size + j]+zt1[i_pos*size + j+1]+zt1[i_pos*size + j-1]-4*zt1[i_pos*size + j]);
			}
		}
	}
}

int main( int argc, char *argv[] ) {
  // verify arguments
  if ( argc != 4 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    return -1;
  }
  int size = atoi( argv[1] );
  int max_time = atoi( argv[2] );
  int interval  = atoi( argv[3] );

  if ( size < 100 || max_time < 3 || interval < 0 ) {
    cerr << "usage: Wave2D size max_time interval" << endl;
    cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
    return -1;
  }
  // create a simulation space
  double z[3][size][size];
  for ( int p = 0; p < 3; p++ ) {
	  for ( int i = 0; i < size; i++ ){
		  for ( int j = 0; j < size; j++ ){
 	      	  z[p][i][j] = 0.0; // no wave
		  }
	  }
  }


  // start a timer
  Timer time;
  time.start( );


  //Initialize pointers
  double (*zt0)[size] =  z[0];
  double (*zt1)[size] =  z[1];
  double (*zt2)[size] =  z[2];

  // time = 0;
  // initialize the simulation space: calculate z[0][][]
  int weight = size / default_size;

  initializeForTimeZero(*zt0, size, weight);

//  cout<<0<<"\n";
//  printCurrentState(*zt0, size);


  // time = 1
  // calculate z[1][][] 
  // cells not on edge
  // IMPLEMENT BY YOURSELF !!!
  //Initialize constants in the equation
  double w1 = (c*c/2)*(dt*dt)/(dd*dd);

  //Loop to initialize all blocks for t=1
  initializeForTimeOne(*zt1, *zt0, size, w1);


//  cout<<1<<"\n";
//  printCurrentState(*zt1, size);


  // simulate wave diffusion from time = 2
  // IMPLEMENT BY YOURSELF !!!

  //Initialize constants in the equation
  double w2 = w1*2;

  //Loop to initialize all blocks for t=2 ....max_time-1
  for ( int t = 2; t <max_time; t++ ){
	  initializeForTimeT(*zt2, *zt1, *zt0, size,w2);
	  //Print current state after gathering values from all ranks
	  if(interval != 0 && (t%interval==0 || t==max_time-1)){
		  cout<<t<<"\n";
		  printCurrentState(*zt2, size);
	  }
	  double (*temp)[size] = zt0;
	  zt0 = zt1;
	  zt1 = zt2;
	  zt2 = temp;
  }
  // end of simulation

  // finish the timer
  cerr << "Elapsed time = " << time.lap( ) << endl;
  return 0;
}
