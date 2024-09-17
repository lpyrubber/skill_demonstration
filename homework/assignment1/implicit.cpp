//exercise 5.3
//method 4, 5
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define L 2.0
#define N 41
#define N_STEP 10
#define N_METHOD 9
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

float *u, *u_new, *v, *gamma, *f;
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
	gamma = (float*)malloc(N_METHOD * N * sizeof(float));
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
		gamma[i] = 0;
		gamma[i + N] = 0;
		v[i] = 0;
		v[i + N] = 0;
	}
	for(i = 1; i < N-1; i++){
		f[i] = u[i];
		f[i + N] = u[i + N] - 0.5 * C * dt / dx * (u[i + N + 1] - u[i + N - 1]);
	}
	//Implicit central
	for(i = 1; i < N - 1; i++){

	}
	A[0] = 1;
	A[N * N - 1] = 1;
	//Crank-Nicolson
	offset = N * N;
	for(i = 1; i < N - 1; i++){
	
	}
}
void Compute_flux(){
	int i, j, offset;
	//boundary value maintain
	bool flag;
	
}
void Update_U(){
	int i, j, offset;
	//calculate new U
	//skip 0 & N for keeping its value
	for(i = 1; i < N-1; i++){
		offset = 0;	//Explicit backward
		u_new[i + offset] = u[i + offset] - C * dt / dx * (u[i + offset] - u[i + offset - 1]);
                offset = N;	//Explicit forward
		u_new[i + offset] = u[i + offset] - C * dt / dx * (u[i + offset + 1] - u[i + offset]);
                offset = 2 * N;	//Explicit central
		u_new[i + offset] = u[i + offset] - 0.5 * C * dt / dx * (u[i + offset + 1] - u[i + offset - 1]);
                offset = 3 * N;	//Lax
		u_new[i + offset] = 0.5 * (u[i + offset + 1] + u[i + offset - 1]) - 0.5 * C * dt / dx * (u[i + offset + 1] - u[i + offset - 1]);
                offset = 4 * N; //Lax-Wendroff
		u_new[i + offset] = u[i + offset] - 0.5 * C * dt / dx * (u[i + offset + 1] - u[i + offset - 1]) + 0.5 * C * C * dt * dt / dx / dx * (u[i + offset + 1] - 2 * u[i + offset] + u[i + offset -1]);
                offset = 5 * N; //MacCormack
		u_new[i + offset] = 0.5 * (u[i + offset] + u_bar[i] - C * dt / dx * (u_bar[i] - u_bar[i - 1]));
                offset = 6 * N; //Jameson
		u_new[i + offset] = u_k[i + 3 * N];
                offset = 7 * N; //Worming-Beam
		if (i == 1){
			u_new[i + offset] = u[i + offset] - C * dt / dx * ((u_h[i] + 0.5 * (u[i + offset] - u[i + offset - 1])) - (u_h[i - 1]));
		}else{	
			u_new[i + offset] = u[i + offset] - C * dt / dx * ((u_h[i] + 0.5 * (u[i + offset] - u[i + offset - 1])) - (u_h[i - 1] + 0.5 * (u[i + offset - 1] - u[i + offset - 2])));
                }
		offset = 8 * N; //Upwind
		u_new[i + offset] = u[i + offset] - dt / dx * (f[i] - f[i - 1]);
	}
	//update to old U
	for(i = 1; i < N-1; i++){
		for(j = 0; j < N_METHOD; j++){
			u[ i + j * N] = u_new[ i + j * N];
		}
	}
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
