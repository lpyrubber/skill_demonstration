#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __APPLE__
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif

#define K 1.0
#define Q 1.0
#define L 1.0
#define N 10
#define dx (L/N)
