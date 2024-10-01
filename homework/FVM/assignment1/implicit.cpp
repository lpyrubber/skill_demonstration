//exercise 5.3
//method 4, 5
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define L 2.0
#define N 41
#define N_STEP 10
#define N_METHOD 2
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

float *u, *u_new, *am, *bm, *cm, *bp, *cp, *dp, *v, *f;
FILE *pFile;

int main(){
	int i, time=0;
	pFile = fopen("implicit_data.txt","w");
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
	f = (float*)malloc(N_METHOD * N * sizeof(float));
	am = (float*)malloc(N_METHOD * N * sizeof(float));
	bm = (float*)malloc(N_METHOD * N * sizeof(float));
	cm = (float*)malloc(N_METHOD * N * sizeof(float));
	bp = (float*)malloc(N * sizeof(float));
	cp = (float*)malloc(N * sizeof(float));
	dp = (float*)malloc(N * sizeof(float));

	u = (float*)malloc(N_METHOD * N * sizeof(float));
	u_new = (float*)malloc(N_METHOD * N * sizeof(float));
}

void Initial(){
	int i = 0, j = 0, offset;
	for( i = 0; i < N; i++){
		if(i*dx < 0.5){
			u[i] = 1;		//Implicit central
			u[i + N] = 1;		//Explicit forward
		}else{
			u[i] = 0.5;		//Implicit backward
			u[i + N] = 0.5;		//Explicit forward
		}
		if(i==0){
			am[i]= 1;
			am[i+N]= 1;
			bm[i]= 0;
			bm[i+N]= 0;
			cm[i]=  0;
			cm[i+N]= 0;
		}else if(i==N-1){
			am[i]= 1 + C * dt / dx;
			am[i+N] = 1 + C * dt / dx;
			bm[i] = 0;
			bm[i+N] = 0;
			cm[i] = -1* C * dt / dx;
			cm[i+N] = -1 * C * dt / dx;
		}else{
			am[i]=1;
			am[i+N]=1;
			bm[i]= 0.5 * C * dt / dx;
			bm[i+N]= 0.25 * C * dt / dx;
			cm[i]= - 0.5 * C * dt / dx;
			cm[i+N]= - 0.25 * C * dt / dx;
		}
	}
	for(i = 1; i < N-1; i++){
		//Implicit central
		f[i] = u[i];
		//Crank-Nicolson
		f[i + N] = u[i + N] - 0.25 * C * dt / dx * (u[i + N + 1] - u[i + N - 1]);
	}
	f[0] = 1;
	f[N] = 1;
	
	f[N - 1] = u[N - 1];
	f[N - 1 + N] = u[N - 1 + N];

}
void Compute_flux(){
	int i, j;

	for(j=0;j<N_METHOD;j++){
		bp[0]=am[0+j*N];
		for(i=1; i<N; i++){
			cp[i]=cm[i+j*N]/bp[i-1];
			bp[i]=am[i+j*N]-cp[i]*bm[i-1+j*N];
		}
		dp[0]=f[0+j*N];
		for(i=1; i<N; i++){
			dp[i]=f[i+j*N]-cp[i]*dp[i-1];
		}
		u_new[N-1+j*N]=dp[N-1]/bp[N-1];
		for(i=N-2; i>=0; i--){
			u_new[i+j*N] = (dp[i]-bm[i+j*N]*u_new[i+1+j*N])/bp[i];
		}
	}
}

void Update_U(){
	int i, j, offset;
	//update old U
	for(i = 0; i < N; i++){
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
	
	f[0] = u[0];
	f[N] = u[N];
	f[N - 1] = u[N-1];
	f[N - 1 + N] = u[N - 1 + N];

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
	free(am);
	free(bm);
	free(cm);
	free(bp);
	free(cp);
	free(dp);
	free(f);
}
