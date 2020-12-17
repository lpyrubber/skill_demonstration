#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>


#define L     1.0
#define N     200
#define dx    (L/N)
#define dt    0.01*dx
#define Z     (dt/dx)
#define gamma 1.4
#define R     1.0
#define CV    (R/(gamma-1))
#define CP    (CV+R)
#define NO_STEP  3200

float *u, *rho, *v, *T;
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
	sprintf(temp, "data.txt");
	pFile = fopen(temp, "w");
	Allocate_Memory();
	Initial();
//	Save_Result();
	for(i=0; i<NO_STEP; i++){
		Compute_Flux();
		Compute_U();
	}
	printf("%f %f\n", u[N], u[N]);
	printf("%f %f %f %f\n", fm[0], fm[1], fp[1], fp[2]);
	Save_Result();
	Free_Memory();
	fclose(pFile);
	return 0;
}

void Allocate_Memory(){
	size_t size=N*sizeof(float);
	rho = (float*)malloc(size);
	v = (float*)malloc(size);
	T = (float*)malloc(size);
	u = (float*)malloc(3*size);
	fm = (float*)malloc(3*size);
	fp = (float*)malloc(3*size);
}

void Initial(){
	int i;
	for(i=0; i<N; i++){
		if(((i-1)<0.5*N)){
			rho[i]=10.0;
			v[i]=0;
			T[i]=1.0;
		}else{
			rho[i]=1.0;
			v[i]=0;
			T[i]=1.0;
		}
		u[i]=rho[i];
		u[i+N]=rho[i]*v[i];
		u[i+N*2]=rho[i]*(CV*T[i]+0.5*v[i]*v[i]);
	}
}

void Compute_Flux(){
	int i; 
	float a;
	float F1, F2,F3;
	float Fr;
	for(i =0; i<N;i++){
		a = sqrt(gamma*R*T[i]);
		Fr = v[i]/a;
		if (Fr >1 ) Fr =1;
		if (Fr < -1) Fr = -1;
		F1 = u[i]*v[i];
		F2 = u[i]*v[i]*v[i] + R*rho[i]*T[i];
		F3 = v[i]*(u[i+N*2] + R*rho[i]*T[i]);
		fp[i] = 0.5*(F1*(Fr+1)+u[i]*a*(1-Fr*Fr));
		fp[i+N] = 0.5*(F2*(Fr+1)+u[i+N]*a*(1-Fr*Fr));
		fp[i+N*2] = 0.5*(F3*(Fr+1)+u[i+N*2]*a*(1-Fr*Fr));
		fm[i] = -0.5*(F1*(Fr-1)+u[i]*a*(1-Fr*Fr));
		fm[i+N] = -0.5*(F2*(Fr-1)+u[i+N]*a*(1-Fr*Fr));
		fm[i+N*2] = -0.5*(F3*(Fr-1)+u[i+N*2]*a*(1-Fr*Fr));
	}
}

void Compute_U(){
	float FL1, FL2, FL3, FR1, FR2, FR3;
	int i;
	for(i=1; i<N-1; i++){
		FL1 = fp[i-1]+fm[i];
		FR1 = fp[i]+fm[i+1];
		FL2 = fp[i-1+N]+fm[i+N];
		FR2 = fp[i+N]+fm[i+1+N];
		FL3 = fp[i-1+N*2]+fm[i+N*2];
		FR3 = fp[i+N*2]+fm[i+1+N*2];
		u[i] = u[i] - Z*(FR1-FL1);
		u[i+N] = u[i+N] - Z*(FR2-FL2);
		u[i+N*2] = u[i+N*2] - Z*(FR3-FL3);
	}
	u[0]   =  u[1];
	u[N]   = -u[1+N];
	u[N*2] =  u[1+N*2];
	u[N-1]     =  u[N-2];
	u[N-1+N]   = -u[N-2+N];
	u[N-1+N*2] =  u[N-2+N*2];
	for(i=0; i<N; ++i){
		rho[i] = u[i];
		v[i] = u[i+N]/u[i];
		T[i] = (u[i+N*2]/u[i]-v[i]*v[i]*0.5)/CV;
	}
}

void Save_Result(){
	int i;
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=1; i<N; i++){
			fprintf(pFile, "%e %e %e %e\n",i*dx, u[i], v[i], T[i]);
		}
		fprintf(pFile,"\n");
	}
}

void Free_Memory(){
	free(fm);
	free(fp);
	free(u);
	free(rho);
	free(v);
	free(T);
}
