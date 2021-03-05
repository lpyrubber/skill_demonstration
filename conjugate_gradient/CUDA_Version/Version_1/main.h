#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#define NO_STEPS 1
#define TPB 32
#define N 64
#define L 1.0
#define dx (L/N)
#define DEBUG 0
#define ERROR 1e-6

void Allocate_Memory(float **h_a, float **h_b, float **h_c);
void Allocate_Device_Memory(float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i, float **d_j);
void Initial(float *h_a, float *h_b, float *h_c);
void Copy_To_Device(float *h_a ,float *h_b, float *h_c, float *d_a, float *d_b, float *d_c, float *d_d, float *d_e, float *d_f, float *d_g);
void Copy_To_Host(float *h_a, float *d_a);
void Save_Result(float *h_a);
void Free_Memory(float **h_a, float **h_b, float **h_c);
void Free_Device_Memory(float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i, float **d_j);
void GPU_VV_Dot(float *d_a, float *d_b, float *d_c, float *d_d, int ele);
void GPU_MV_Dot(float *d_a, float *d_b, float *d_c, int ele);
void GPU_Update_P(float *d_a, float *d_b, float *d_c, float *d_d, int ele);
void GPU_Update_S(float *d_a, float *d_b, float *d_c, float *d_d, int ele);
void GPU_Update_X(float *d_a, float *d_b, float *d_c, float *d_d, int ele);
void GPU_Update_R(float *d_a, float *d_b, float *d_c, float *d_d, int ele);
void GPU_Update_temp(float *h_a, float *d_a);
