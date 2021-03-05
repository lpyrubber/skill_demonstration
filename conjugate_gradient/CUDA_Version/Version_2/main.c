#include "main.h"
#include <math.h>

int main(){
	float *A, *b, *x, temp;
	float *d_A, *d_b, *d_x, *d_p, *d_s, *d_r, *d_temp, *d_scalar;
	int i;
	char flag;
	Allocate_Memory( &A, &b, &x, &d_A, &d_b, &d_x \
		       , &d_p, &d_s, &d_r, &d_temp, &d_scalar);
	Initial(A, b, x);
	Send_To_Device(A, b, x , d_A, d_b, d_x);
	i = 0;
	flag = 1;
	Device_Initial(d_A, d_b, d_x, d_r, d_p, d_temp, d_scalar);
	while( ( flag )&&( i < NO_STEP ) ){
		Conjugate_Gradient( d_A, d_b, d_x, d_p, d_s \
				  , d_r, d_temp, d_scalar, &temp );
		if(sqrt(temp) < ERROR ){
			flag=0;
			printf("break at %d\n", i);
		}
		printf("Error = %e\n", temp);
		i++;
	}
	Send_To_Host(x, d_x);
	Save_Result(x);
	Free_Memory( &A, &b, &x, &d_A, &d_b, &d_x \
		   , &d_p, &d_s, &d_r, &d_temp, &d_scalar);
	return 0;
}
