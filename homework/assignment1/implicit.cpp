#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 2.0
#define N 41
#define N_STEP 10
#define N_METHOD 2 
#define dx L/(N-1)
#define C 1
#define CFL 0.9
#define dt CFL*dx/C
