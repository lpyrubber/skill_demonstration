//exercise 5.1
//method 1, 2, 3, 6, 7, 8,  9, 10, 11
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define L 2.0
#define N 41
#define N_STEP 10
#define N_METHOD 9
#define dx (float)(L/(N-1))
#define C 1
#define CFL 1.0
#define dt (float)(CFL*dx/C)

void Allocate_Memory();
void Initial();
void Compute_flux();
void Update_U();
void Save_Result(int time);
void Free_Memory();

float *u, *u_new, *u_bar, *u_k, *u_h, *f;
FILE *pFile;

int main(){
	int i, time=0;
	pFile = fopen("explicit_data.txt","w");
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
	u = (float*)malloc(N_METHOD * N * sizeof(float));
	u_new = (float*)malloc(N_METHOD * N * sizeof(float));
	u_bar = (float*)malloc(N * sizeof(float));
	u_k = (float*)malloc(4 * N * sizeof(float));
	u_h = (float*)malloc(N * sizeof(float));
	f = (float*)malloc(N * sizeof(float));
}

void Initial(){
	int i = 0;
	for( i = 0; i < N; i++){
		if(i*dx < 0.5){
			u[i] = 1;		//Explicit backward
			u[i + N] = 1;		//Explicit forward
			u[i + 2 * N] = 1;	//Explicit central
			u[i + 3 * N] = 1;	//Lax
			u[i + 4 * N] = 1;	//Lax-Wendroff
			u[i + 5 * N] = 1;	//MacCormack
			u[i + 6 * N] = 1;	//Jameson
			u[i + 7 * N] = 1;	//Worming-Beam
			u[i + 8 * N] = 1;	//Upwind
		}else{
			u[i] = 0.5;		//Explicit backward
			u[i + N] = 0.5;		//Explicit forward
			u[i + 2 * N] = 0.5;	//Explicit central
			u[i + 3 * N] = 0.5;	//Lax
			u[i + 4 * N] = 0.5;	//Lax-Wendroff
			u[i + 5 * N] = 0.5;	//MacCormack
			u[i + 6 * N] = 0.5;	//Jameson
			u[i + 7 * N] = 0.5;	//Worming-Beam
			u[i + 8 * N] = 0.5;	//Upwind
		}
	}
}
void Compute_flux(){
	int i, j, offset;
	//boundary value maintain
	u_bar[0] = 1;
	u_bar[N - 1] = 0.5;
	for( j = 0; j < 4; j++){
		u_k[j * N] = 1;
		u_k[N - 1 + j * N] = 0.5; 	
	}
	u_h[0] = 1;
	u_h[N - 1] = 0.5;
	f[0] = C * 1;
	f[N - 1] = C * 0.5;
	for( j = 0; j < 4; j ++){
		for( i = 1; i < N-1; i++){
			if(j == 0){
				//MacCormack
				offset = 5 * N;
				u_bar[i] = u[i + offset] - C * dt / dx * (u[i + offset + 1] - u[i + offset]);  
				//Jameson
				offset = 6 * N;
				u_k[i] = u[i + offset] - 0.125 * C * dt / dx * (u[i + offset + 1] - u[i + offset - 1]);
				//warming-Beam
				offset = 7 * N;
				u_h[i] = u[i + offset] - 0.5 * C * dt / dx * (u[i + offset] - u[i + offset - 1]);
				//Upwind
				offset = 8 * N;
				f[i] = 0.5 * C * (u[i + offset + 1] + u[i + offset]) - 0.5 * fabs(C) * (u[ i + offset + 1] - u[i + offset]);
			}else{
				//Jameson
				offset = 6 * N;
				u_k[i + j * N] = u[i + offset] - 0.5 * C * dt / dx / (4 - j) * (u_k[i + (j - 1) * N + 1] - u_k[i + (j - 1) * N - 1]);
			}
		}
	}
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
