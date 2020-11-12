//load header file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


//define constant value
#define L       100.0       //miles of total length 
#define N       200        //total number of dx
#define DX      (L/N)      //miles of dx
#define DT      (0.01*DX)  //sec for dt
#define Z       (DT/DX)    //
#define GAP     10	   //pick the value after number of dx
#define NO_STEP 800         //total step
#define G       9.81       //gravity acceleration
#define DEBUG   1

//define constant value for MPI
#define NONE   0
#define MASTER 0
#define LEFT   1
#define RIGHT  0
#define BEGIN  0
#define FINAL  1
//define constant for openMP
#define NT     2

//preclaim function
void Allocate_Memory( int id , int max , float **a , float **b , float **c );
void Compute( int id , int ele , float *a , float *b , float *c );
void Free_Memory( float **a , float **b , float **c );
