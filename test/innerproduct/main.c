#include "main.h"

int main(){
	float *A, *dot;
	float *d_A, *d_dot, sum=0;
	int i, num;
	Allocate_Memory( &A, &dot , &d_A, &d_dot);
	Initial( A );
	Send_To_Device( A , d_A );
	GPU_Dot( d_A , d_dot );
	Send_To_Host( dot, d_dot );
	for( i = 0; i < BPG ; i++){
		sum += dot[i];
		num = ( i == BPG - 1 ) ? N-(BPG-1)*TPB : TPB;
		printf("%d: %d  %e %e\n", i , num, dot[i], dot[i]/num);
	}
	printf("result = %f\n", dot[0]/N);
	Free_Memory( &A, &dot, &d_A, &d_dot );
	return 0;
}

