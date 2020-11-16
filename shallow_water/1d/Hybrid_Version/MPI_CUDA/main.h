//load header file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//define constant value
#define L       100.0       //miles of total length 
#define N       200        //total number of dx
#define DX      (L/N)      //miles of dx
#define DT      (0.01*DX)  //sec for dt
#define Z       (DT/DX)    //
#define GAP     10	   //pick the value after number of dx
#define NO_STEP 800         //total step
#define G       9.81       //gravity acceleration
#define DEBUG   0

//define constant value for MPI
#define NONE   0
#define MASTER 0
#define LEFT   1
#define RIGHT  0
#define BEGIN  0
#define FINAL  1

//define constant value for cuda
#define TPB     256
#define BPG     (int)((N+TPB-1)/TPB )

//preclaim function
void Allocate_Memory( int id , int max , float **h_a , float **h_b , float **d_a , float **d_b , float **d_c );
void Sent_To_Device( int ele , float *h_a , float *h_b , float *d_a , float *d_b );
void Compute( int ele , float *d_a , float *d_b , float *d_c );
void Sent_To_Host( int ele , float *h_a , float *h_b , float *d_a , float *d_b );
void Free_Memory( int id , float **h_a , float **h_b , float **d_a , float **d_b , float **d_c );
