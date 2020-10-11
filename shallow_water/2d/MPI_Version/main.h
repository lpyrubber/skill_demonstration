//load header
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//define constant value
#define Lx      100.0
#define Ly      100.0
#define Nx      200
#define Ny	200
#define DX	(Lx/Nx)
#define DY      (Ly/Ny)
#define DT      (0.01*DX)
#define Zx      (DT/DX)
#define Zy      (DT/DY)
#define NO_STEP 800
#define G       9.81
#define DEBUG   1

//define constant for MPI
#define NONE   0
#define MASTER 0
#define LEFT   1
#define RIGHT  0
#define BEGIN  0
#define FINAL  1

//preclaim function
void Allocate_Memory( int id , int max , float **a , float **b , float **c );
void Compute( int id , int ele , float *a , float *b , float *c );
void Free_Memory( float **a , float **b , float **c );
