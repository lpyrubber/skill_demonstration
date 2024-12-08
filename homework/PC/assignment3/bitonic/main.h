#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define TPB 256


void Allocate_Memory();
void Send_To_Device();
void GPU_Compute();
void Send_To_Host();
void Free_Memory();


extern int *h_a, *d_a;
extern int N, N2;
extern int BPG; 