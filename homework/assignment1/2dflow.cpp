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
#define N_it 20
#define NO (int)(0.5*(NX-NC))
#define L (NX-1)*C
#define D (NY-1)*C
#define LO (L-C)
#define dxmin 0.05*C
#define dymin 0.1*TH
#define N_METHOD 5

void Allocate_Memory();
void Initial();
void Create_Grid();
//void Interation();
float Calculate_Kappa(int N, float length, float delta);
void Save_Result(int time);
void Free_Memory();

float *x, *y, *A, *phi, *phi_new, *b;
FILE *pFile;

int main(){
	int i, time=0;
	pFile = fopen("data.txt","w");
	Allocate_Memory();
	Create_Grid();
	Initial();
	Save_Result(time);
//	Interation();
	time++;
//	Save_Result(time);
	fclose(pFile);
	Free_Memory();

}


void Allocate_Memory(){
	x = (float*)malloc(NX * sizeof(float));
	y = (float*)malloc(NY * sizeof(float));
//	A = (float*)malloc(NX * NX * NY * NY * sizeof(float));
	phi = (float*)malloc(NX * NY * N_METHOD * sizeof(float));
	phi_new = (float*)malloc(NX * NY * N_METHOD * sizeof(float));
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
	for(i = 0; i < NX; i++)
		printf("%2.2f ",x[i]);
	K= Calculate_Kappa(NY, D, dymin);
	y[0]=0;
	for(i = 1; i < NY; i++){
		y[i]=D*(exp(K*(i)/(NY-1))-1)/(exp(K)-1);
	}

}

float Calculate_Kappa(int N, float length, float delta){
	float K=1, i;
	for(i=0; i<N_it; i++){
		K=K-(delta-length*(exp(K/(N-1))-1)/(exp(K)-1))/(length*((N-1)*(exp(K/(N-1))-1)*exp(K)-(exp(K)-1)*exp(K/(N-1)))/((N-1)*(exp(K)-1)*(exp(K)-1)));
	}
	printf("K=%e\n",K);
	return K;
}

void Initial(){
	int i, j, k;
	for(k=0; k<N_METHOD; k++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				phi[i+j*NX+k*NX*NY]=0;
			}
		}
	}
	//set up background
	for(i=0; i<NC; i++){
		phi[NO+i+1]=phi[NO+i-1]+(x[NO+i+1]-x[NO+i-1])*(-(x[NO+i]-0.5*C)/sqrt(pow((0.25*(C*C)/TH+0.25*TH),2)-pow((x[NO+i]-0.5*C),2)));	
	}
}
void Save_Result(int time){
	int i, j, k;
	if(time == 0){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				fprintf(pFile, "%e ", x[i]);
			}
		}
		fprintf(pFile, "\n");
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				fprintf(pFile, "%e ", y[i]);
			}
		}
		fprintf(pFile, "\n");
		for(j = 0; j < NY; j++){ 
			for(i = 0; i < NX; i++){
				fprintf(pFile, "%e ", phi[i + j * NX]);	
			}
		}
		fprintf(pFile, "\n");
	}else{
		for(k = 0; k < N_METHOD; k++){
			for(j = 0; j < NY; j++){
				for( i = 0; i < NX; i++){
					fprintf(pFile, "%e ", phi[ i + j * NX + k * NX * NY ]);
				}
				fprintf(pFile, "\n");
			}
		}
	}
}

void Free_Memory(){
	free(phi);
	free(phi_new);
	free(x);
	free(y);
	free(b);
//	free(A);
}
