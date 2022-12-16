#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define L       100.0
#define N       20000
#define DX      (L/N)
#define DT      (0.01*DX)
#define Z       (DT/DX)
#define TPB     256
#define BPG     (int)((N+TPB-1)/TPB )
#define NO_STEP 80000
#define G       9.81

#define DEBUG 0

void Allocate_Memory( float **h_a , float **h_b , float **d_a , float **d_b , float **d_c );
void Sent_To_Device( float *h_a , float *h_b , float *d_a , float *d_b );
void GPU_Compute( float *d_a , float *d_b , float *d_c );
void Sent_To_Host( float *h_a , float *h_b , float *d_a , float *d_b );
void Free( float **h_a , float **h_b , float **d_a , float **d_b , float **d_c );
