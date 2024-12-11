#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define N_IT 20
#define SUM_MAX 1e14
#define USE_MATRIX 1

//c function
void Allocate_Memory();
void Free_Memory();
void Send_To_Device();
void GPU_Compute();
void Send_To_Host();
void Save_Result();

extern int NC, NC2, TPB, BPG, BPG2, BPGC, NP, NP2, Dim, NT, NB;

extern int *h_label, *h_id, *h_id_old, *h_flag;
extern float *h_x;

extern int *d_label, *d_id, *d_id_old, *d_nlist, *d_prefix, *d_flag;
extern float *d_x, *d_min, *d_sum;