#include <stdlib.h>
#include <stdio.h>

#define DEBUG 0
#define K 1.0
#define Q 1.0
#define L 1.0
#define N 64
#define dx (L/N)
#define NO_STEP 2
#define ERROR 1e-6
#define TPB 4
#define BPG ((N+TPB-1)/TPB)

void Allocate_Memory( float **h_a, float **h_b, float **h_c, float **d_a, float **d_b, float **d_c \
		    , float **d_d, float **d_e, float **d_f, float **d_g, float **d_h );
void Initial( float *a, float *b, float *c );
void Send_To_Device( float *h_a, float *h_b, float *h_c, float *d_a, float *d_b , float *d_c );
void Device_Initial( float *a, float *b, float *c, float *d, float *e, float *f, float *g );
void Conjugate_Gradient( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h ,float *temp);
void Send_To_Host( float *h_a, float *d_a );
void Save_Result( float *a );
void Free_Memory( float **h_a, float **h_b, float **h_c, float **d_a, float **d_b, float **d_c \
		, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h );
