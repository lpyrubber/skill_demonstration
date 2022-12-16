#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define L 100.0
#define N 20000
#define DX (L/N)
#define DT (0.01*DX)
#define Z (DT/DX)
#define GAP 10
#define NO_STEP 80000
#define G 9.81

float *u, *u_new;
float *fm, *fp;
FILE *pFile;
void Allocate_Memory();
void Initial();
void Compute_Flux();
void Compute_U();
void Update_U();
void Save_Result();
void Free_Memory();

int main(){
	static int k=0;
	char temp[20];
	int i;
	double total_t;
	clock_t start_t, end_t;
	start_t=clock();
	sprintf(temp, "data_%d.txt",k);
	pFile = fopen(temp, "w");
	Allocate_Memory();
	Initial();
//	Save_Result();
	for(i=0; i<NO_STEP; i++){
		Compute_Flux();
		Compute_U();
		Update_U();
//		if(i%GAP==(GAP-1)&&i<NO_STEP-1){
//			Save_Result();
//		}
	}
//	printf("%f %f\n", u[N+2], u[1+N+2]);
//	printf("%f %f %f %f\n", fm[0], fm[1], fp[1], fp[2]);
	Save_Result();
	Free_Memory();
	fclose(pFile);
	end_t=clock();
	total_t=(double)(end_t-start_t)/CLOCKS_PER_SEC;
	printf("Total time taken by CPU: %lf sec\n", total_t);
	return 0;
}

void Allocate_Memory(){
	size_t size=2*(N+2)*sizeof(float);
	u = (float*)malloc(size);
	u_new = (float*)malloc(size);
	fm = (float*)malloc(size);
	fp = (float*)malloc(size);
}

void Initial(){
	int i;
	for(i=0; i<N+2; i++){
		if(((i-1)<0.5*N)){
			u[i]=10.0;
			u[i+N+2]=0;
		}else{
			u[i]=1.0;
			u[i+N+2]=0;
		}
		
//		u[i]=10.0;
//		u[i+N+2]=0;
	}
}

void Compute_Flux(){
	int i; 
	float vel,a;
	float F1, F2;
	float Fr;
	for(i =0; i<N+2;i++){
		vel=u[i+N+2]/u[i];
		a = sqrt(G*u[i]);
		Fr = vel/a;
		if (Fr >1 ) Fr =1;
		if (Fr < -1) Fr = -1;
		F1 = u[i]*vel;
		F2 = u[i]*vel*vel + 0.5*G*u[i]*u[i];
		fp[i] = 0.5*(F1*(Fr+1)+u[i]*a*(1-Fr*Fr));
		fp[i+N+2] = 0.5*(F2*(Fr+1)+u[i+N+2]*a*(1-Fr*Fr));
		fm[i] = -0.5*(F1*(Fr-1)+u[i]*a*(1-Fr*Fr));
		fm[i+N+2] = -0.5*(F2*(Fr-1)+u[i+N+2]*a*(1-Fr*Fr));
	}
}

void Compute_U(){
	float FL1, FL2, FR1, FR2;
	int i;
	for(i=1; i<N+1; i++){
		FL1 = fp[i-1]+fm[i];
		FR1 = fp[i]+fm[i+1];
		FL2 = fp[i+N-1+2]+fm[i+N+2];
		FR2 = fp[i+N+2]+fm[i+N+1+2];
		u_new[i] = u[i] - Z*(FR1-FL1);
		u_new[i+N+2] = u[i+N+2] - Z*(FR2-FL2);
	}
	u_new[0] = u_new[1];
	u_new[N+2] = -u_new[1+N+2];
	u_new[N+1] = u_new[N];
	u_new[N+1+N+2] = -u_new[N+N+2];
}

void Update_U(){
	int i;
	for(i=0; i<N+2; i++){
		u[i] = u_new[i];
		u[i+N+2] = u_new[i+N+2];
	}
}

void Save_Result(){
	int i;
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=1; i<N+1; i++){
			fprintf(pFile, "%e %e %e\n",i*DX, u[i], u[i+N+2]);
		}
		fprintf(pFile,"\n");
	}
}

void Free_Memory(){
	free(fm);
	free(fp);
	free(u);
	free(u_new);
}
