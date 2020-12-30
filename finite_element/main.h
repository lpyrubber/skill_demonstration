#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NU    0.3    //Poisson's ratio
#define E     210e9   //Modulus of Elasticity
#define THICK 0.025  //Thickness of plate

void Count_Number(int *np, int*ne );
void Allocate_Memory(int np, int ne, char **Fix_x, char **Fix_y, int **tri \
		    , float **x, float **y, float **Fx, float **Fy, float **K \
		    , float **Kr, float **Fr, float **u, float **U, float **S );
void Initial(int np, int ne, char *Fix_x, char *Fix_y, int *tri, float *x, float *y, float *Fx, float *Fy);
void Calculate_K(int np, int ne, int *tri, float *x, float *y , float *K);
void Partition_K_F( int np, int *nr, char *Fix_x, char *Fix_y, float *Fx, float *Fy \
		   , float *Fr, float *K, float *Kr);
void Calculate_Stress();
void Partition_U(int np, int nr, char *Fix_x , char *Fix_y , float *u, float *U);
void Free_Memory( char **Fix_x, char **Fix_y, int **tri \
		, float **x, float **y, float **Fx, float **Fy, float **K \
		, float **Kr, float **Fr, float **u, float **U, float **S );
void Conjugate_Gradient(float *A, float *b, float *x, int N);
void Save_Result(int np, int ne, int *tri, float *a, float *b);
