#include "main.h"

int main(){
	float *h_x;
	float *d_x, *d_scalar;
	int i;
	Allocate_Memory(&h_x, &d_x, &d_scalar);
	Compute(d_x, d_scalar);
	Send_To_Host(h_x, d_x);
	for( i = 0 ; i < N ; i++){
		printf("%f\n", h_x[i]);
	}
	Free_Memory(&h_x, &d_x, &d_scalar );
	return 0;
}
