#include <stdio.h>
#include <stdlib.h>

#define TPB     256
#define BPG     (int)(((Nx+2)*(Ny+2)+TPB-1)/TPB )
#define Lx      100.0
#define Ly      100.0
#define Nx      200
#define Ny      200
#define	DX      (Lx/Nx)
#define DY      (Ly/Ny)
#define DT      (0.01*DX)
#define Zx      (DT/DX)
#define Zy      (DT/DY)
#define NO_STEP 800
#define G       9.81

#define DEBUG   0

void Allocate_Memory( float **h_a , float **d_a , float **d_b , float **d_c );
void Sent_To_Device( float *h_a , float *d_a );
void GPU_Compute( float *d_a , float *d_b , float *d_c );
void Sent_To_Host( float *h_a , float *d_a );
void Free( float **h_a , float **d_a , float **d_b , float **d_c );
