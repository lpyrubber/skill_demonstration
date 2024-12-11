#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_points 10

int main(){
	int i, ix,iy;
	int N=0.5*(N_points+1)*N_points;
	int count=N_points;
	int k=0;
	int add=N_points-1;
	double temp;
/*	printf("N=%d\n",N);
	for(i=0;i<N;i++){
		ix=0.5*((2*N_points-1)-sqrtf(4*N_points*N_points+4*N_points+1-8*(i)))+1;
		temp=0.5*((2*N_points-1)-sqrtf(4*N_points*N_points+4*N_points+1-8*(i)))+1;
		iy=i-0.5*(2*N_points-ix+1)*(ix);
		printf("%d, %d, temp=%f\n",ix,ix+iy, temp);
		if(i==add){
			count--;
			add+=count;
			printf("next line\n");
		}
	}
*/	
	int id;
	printf("N=%d\n",N_points);
	for(i=0;i<N_points;i++){
		id=(i%2)?((i%2)*(N_points-1)-i/2):((i%2)*(N_points-1)+i/2);
		printf("i=%d, ib=%d\n",i, id);
	}
	printf("\n");
	return 0;
}
