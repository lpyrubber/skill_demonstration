#include <stdlib.h>
#include <stdio.h>

#define N 100000
#define TPB 256
#define BPG ( ( int )( N + TPB - 1 ) / TPB )
#define DEBUG 0


void Allocate_Memory(float **h_a, float **h_b, float **d_a, float **d_b);
void Initial(float *h_a );
void Send_To_Device( float *h_a , float *d_a );
void GPU_Dot( float *a, float *b );
void Send_To_Host( float *h_a, float *d_a );
void Free_Memory( float **h_a, float **h_b, float **d_a, float **d_b );

