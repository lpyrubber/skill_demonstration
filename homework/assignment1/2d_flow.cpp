#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C 1.0
#define TH 0.06
#define NY 51
#define NX 51
#define L (NX-1)*C
#define D (NY-1)*C
#define dy_min 0.1*TH
#define dx_min 0.05*C
#define NX_out 0.5*(NX-21)+1
