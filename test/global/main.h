#include <stdlib.h>
#include <stdio.h>

#define N 64
#define TPB 32
#define BPG ( ( N + TPB - 1 ) / TPB )

void Allocate_Memory( float **h_a, float **d_a, float **d_b );
void Compute( float *a, float *b );
void Send_To_Host( float *h_a, float *d_a );
void Free_Memory( float **h_a, float **d_a, float **d_b );
