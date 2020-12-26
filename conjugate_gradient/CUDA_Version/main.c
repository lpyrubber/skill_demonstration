#include "main.h"

int main(){
	float *h_A, *h_x, *h_b;
	float *d_A, *d_x, *d_b, *d_p, *d_s, *d_t, *d_r, *d_v, *d_temp, *d_scalar;
	float temp;
	int i;	
	Allocate_Memory(&h_A, &h_x, &h_b);
	Allocate_Device_Memory(&d_A, &d_x, &d_b, &d_p, &d_s, &d_t, &d_r, &d_v, &d_temp, &d_scalar);
	Initial(h_A, h_x, h_b);
	Copy_To_Device(h_A, h_x, h_b, d_A, d_x, d_p, d_v, d_b, d_r, d_scalar);
	//start compute
	for(int i=0; i<NO_STEPS; i++){
		GPU_VV_Dot(d_r, d_b, d_temp, d_scalar, N);
		GPU_Update_P(d_p, d_r, d_v, d_scalar, N);
		GPU_MV_Dot(d_A, d_p, d_v, N);
		GPU_VV_Dot(d_b, d_v, d_temp, d_scalar+6, N);
		GPU_Update_S(d_s, d_r, d_v, d_scalar, N);
		GPU_MV_Dot(d_A, d_s, d_t, N);
		GPU_VV_Dot(d_t, d_s, d_temp, d_scalar+6, N);
		GPU_VV_Dot(d_t, d_t, d_temp, d_scalar+7, N);
		GPU_Update_X(d_x, d_p, d_s, d_scalar, N);
		GPU_Update_R(d_r, d_s, d_t, d_scalar, N);
		GPU_VV_Dot(d_r, d_r, d_temp, d_scalar+6, N);
		GPU_Update_temp(&temp, d_scalar+6);
		if(sqrt(temp)<ERROR){
			printf("break at %d\n",i);
			break;
		}
	}
	Copy_To_Host(h_x, d_x);
	Save_Result(h_x);
	Free_Memory(&h_A, &h_x, &h_b);
	Free_Device_Memory(&d_A, &d_x, &d_b, &d_p, &d_s, &d_t, &d_r, &d_v, &d_temp, &d_scalar);
}
