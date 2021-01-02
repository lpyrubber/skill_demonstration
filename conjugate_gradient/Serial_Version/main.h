#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define K 1.0
#define Q 1.0
#define L 1.0
#define N 64
#define dx (L/N)

void Conjugate_Gradient(float *A, float *b, float *x);
void Allocate_Memory( float **A, float **b, float **x);
void Initial( float *A, float *b, float *x);
void Save_Result( float *x );
void Free_Memory( float **A, float **b, float **x );
