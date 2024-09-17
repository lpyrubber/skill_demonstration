//exercise 5.3
//method 4, 5
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define L 2.0
#define N 41
#define N_STEP 10
#define N_METHOD 2
#define dx L/(N-1)
#define C 1
#define CFL 0.9
#define dt CFL*dx/C

void Allocate_Memory();
void Initial();
void Compute_flux();
void Update_U();
void Save_Result(int time);
void Free_Memory();

float *u, *u_new, *Gamma, *v, *f;
FILE *pFile;

int main(){
	int i, time=0;
	pFile = fopen("data.txt","w");
	Allocate_Memory();
	Initial();
	Save_Result(time);
	for(i = 0; i < N_STEP; i++){
		Compute_flux();
		Update_U();
	}
	time++;
	Save_Result(time);
	fclose(pFile);
	Free_Memory();

}


void Allocate_Memory(){
	v = (float*)malloc(N_METHOD * N * sizeof(float));
	Gamma = (float*)malloc(N_METHOD * N * sizeof(float));
	f = (float*)malloc(N_METHOD * N * sizeof(float));
	u = (float*)malloc(N_METHOD * N * sizeof(float));
	u_new = (float*)malloc(N_METHOD * N * sizeof(float));
}

void Initial(){
	int i = 0, j = 0, offset;
	for( i = 0; i < N; i++){
		if(i*dx < 0.5){
			u[i] = 1;		//Explicit backward
			u[i + N] = 1;		//Explicit forward
		}else{
			u[i] = 0.5;		//Explicit backward
			u[i + N] = 0.5;		//Explicit forward
		}
		Gamma[i] = 0;
		Gamma[i + N] = 0;
		v[i] = 0;
		v[i + N] = 0;
	}
	for(i = 1; i < N-1; i++){
		//Implicit central
		f[i] = u[i];
		//Crank-Nicolson
		f[i + N] = u[i + N] - 0.25 * C * dt / dx * (u[i + N + 1] - u[i + N - 1]);
	}
	
	f[1] -= -0.5 * C * dt / dx * u[0];
	f[N + 1 ] -= -0.25 * C * dt /dx * u[N];

}
void Compute_flux(){
	int i, j;

	Gamma[N - 1] = 0; 
	Gamma[N + N - 1] = 0;
	for (i = N - 2; i > 0; i--){
		Gamma[i] = - 0.5 * C * dt / dx / ( 1 - 0.5 * C * dt * Gamma[i]);
		Gamma[i + N] = - 0.25 * C * dt / dx / ( 1 - 0.25 * C * dt * Gamma[i + N + 1]);
		
		v[i] = (f[i] - 0.5 * C * dt / dx * v[i + 1]) / (1 - 0.5 * C * dt * Gamma[i + 1]) ;
		v[i + N] = (f[i + N] - 0.25 * C * dt / dx * v[i + N + 1]) / (1 - 0.25 * C * dt * Gamma[i + N + 1]) ;
	}
	u_new[1] = v[1];
	for (i = 2; i < N - 1; i++){
		u_new[i] = v[i] - Gamma[i] * u_new[i - 1];
		u_new[i + N] = v[i + N] - Gamma[i + N] * u_new[i + N - 1];
	}
}
void Update_U(){
	int i, j, offset;
	//update old U
	for(i = 1; i < N-1; i++){
		for(j = 0; j < N_METHOD; j++){
			u[ i + j * N] = u_new[ i + j * N];
		}
	}
	//update f
	for(i = 1; i < N-1; i++){
		//Implicit central
		f[i] = u[i];
		//Crank-Nicolson
		f[i + N] = u[i + N] - 0.25 * C * dt / dx * (u[i + N + 1] - u[i + N - 1]);
	}
	
	f[1] -= -0.5 * C * dt / dx * u[0];
	f[N + 1 ] -= -0.25 * C * dt /dx * u[N];

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
	free(Gamma);
	free(v);
	free(f);
}
