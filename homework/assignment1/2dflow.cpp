//exercise 5.1
//method 1, 2, 3, 6, 7, 8,  9, 10, 11
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Minf 0.5
#define C 1.0
#define TH 0.06
#define Pinf 1.0
#define Rhoinf 1.0 
#define NX 51
#define NY 51
#define NC 21
#define NO 0.5*(NX-NC)
#define L (NX-1)*C
#define D (NY-1)*C
#define LO L-C
#define dxmin 0.05*C
#define dymin 0.1*TH
#define N_METHOD 5

void Allocate_Memory();
void Initial();
void Create_Grid();
void Interation();
void Calculate_Kappa(int N, float length, float delta);
void Save_Result(int time);
void Free_Memory();

float *x, *y, *A, *phi, *phi_new, *f;
FILE *pFile;

int main(){
	int i, time=0;
	pFile = fopen("data.txt","w");
	Allocate_Memory();
	Create_Grid();
	Initial();
	Save_Result(time);
	Interation();
	time++;
	Save_Result(time);
	fclose(pFile);
	Free_Memory();

}


void Allocate_Memory(){
	x = (float*)malloc(NX * sizeof(float);
	y = (float*)malloc(NY * sizeof(float);
	A = (float*)malloc(NX * NX * NY * NY * sizeof(float));
	phi = (float)malloc(NX * NY * N_METHOD * sizeof(float));
	phi_new = (float)malloc(NX * NY * N_METHOD * sizeof(float));
	b = (float*)malloc(NX * NY * N_METHOD * sizeof(float));
}

void Create_Grid(){
	int i;
	float K;
	for( i = 0; i < NC; i++){
		x[NO + i] = i * dxmin;
	}
	K = Calculate_Kappa(NO+1, LO, dxmin);
	for(i = 0; i < NO; i++){
		x[NO-1-i] = -LO*(exp(K*(i+1)/NO)-1)/(exp(K)-1);
		x[NO + NC + i] = C + LO*(exp(K*(i+1)/NO)-1)/(exp(K)-1); 
	}
	K= Calculate_Kappa(NY, D, dymin);
	for(i = 0; i < NY; i++){
		y[i]

}

void Save_Result(int time){
	int i, j;
	if(time == 0){
		for(i = 0; i < N; i++){
			fprintf(pFile, "%e ", i * dx);
		}
		fprintf(pFile, "\n");
		for(i = 0; i < N; i++){
			fprintf(pFile, "%e ", u[i]);	
		}
		fprintf(pFile, "\n");
	}else{
		for(j = 0; j < N_METHOD; j++){
			for( i = 0; i < N; i++){
				fprintf(pFile, "%e ", u[ i + j * N]);
			}
			fprintf(pFile, "\n");
		}
	}
}

void Free_Memory(){
	free(u);
	free(u_new);
	free(u_bar);
	free(u_k);
	free(u_h);
	free(f);
}
*/
