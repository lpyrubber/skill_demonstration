#include <stdio.h>
#include <stdlib.h>

#define L 1.0
#define N 200
#define DX (L/N)
#define DT (0.01*DX)
#define Z (DT/DX)
#define NO_STEP 3200
#define R 1.0
#define gamma 1.4
#define CV (R/(gamma-1))
#define CP (CV+R)

#define TPB 256
#define BPG ((int)((N+TPB-1)/TPB))
#define DEBUG 0

void Allocate_Memory( float **h_a, float **h_b, float **h_c, float **h_d\
		    , float **d_a, float **d_a, float **d_c, float **d_d, float **d_e, float **d_f );
void Send_To_Device( float *h_a, float *h_b, float *h_c, float *h_d\
		   , float *d_a, float *d_b, float *d_c, float *d_d);
void GPU_Compute( float *d_a, float *d_b, float *d_c, float *d_d, float *d_e, float *d_f );
void Send_To_Host( float *h_a, float *h_b, float *h_c, float *h_d\
		 , float *d_a, float *d_b, float *d_c, float *d_d);
void Free( float **h_a, **h_b, **h_c, **h_d\
	 , float **d_a, **d_b, **d_c, **d_d, **d_e, **d_f );
