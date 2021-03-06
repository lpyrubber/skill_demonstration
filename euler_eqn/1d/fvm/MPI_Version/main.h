//load header file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//define constant value
#define L       1.0       //miles of total length 
#define N       200        //total number of dx
#define DX      (L/N)      //miles of dx
#define DT      (0.01*DX)  //sec for dt
#define Z       (DT/DX)    //
#define gamma 1.4
#define R     1.0
#define CV    (R/(gamma-1))
#define CP    (CV+R)
#define NO_STEP  3200
#define DEBUG   1

//define constant value for MPI
#define NONE   -1
#define MASTER 0
#define LEFT   1
#define RIGHT  0
#define BEGIN  0
#define FINAL  1

//preclaim function
void Allocate_Memory( int id , int max , float **a , float **b , float **c , float **d , float **e , float **f);
void Compute( int id , int ele , float *a , float *b , float *c , float *d , float *e , float *f);
void Free_Memory( float **a , float **b , float **c , float **d , float **e , float **f);
