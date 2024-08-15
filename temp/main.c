#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>


#define N_SPECIES 1
#define L     1.0
#define N     200
#define Tdrop 20.0
#define Tamp  100.0
#define dx    (L/N)
#define dt    0.01*dx
#define Z     (dt/dx)
#define gamma 1.4
#define R     1.0
#define CV    (R/(gamma-1))
#define CP    (CV+R)
#define NO_STEP  3200
#define mu 1.0
#define lamda 1.0

float *u, *rho, *v, *T, *D, *Y, *x, *W;
float *fm, *fp;
FILE *pFile;
void Allocate_Memory();
void Initial();
void Compute_Flux();
void Compute_U();
void Update_U();
void Save_Result();
void Free_Memory();
void Compute_D();

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
	size_t size=N_SPECIES*N*sizeof(float);

	rho = (float*)malloc(size);
	v = (float*)malloc(size);
	Y = (float*)malloc(size);
	x = (float*)malloc(size);
	T = (float*)malloc(size);
	u = (float*)malloc(3*size);
	fm = (float*)malloc(3*size);
	fp = (float*)malloc(3*size);
	D = (float*)malloc(N_SPECIES*sizeof(float));
	W = (float*)malloc(N_SPECIES*sizeof(float));
}

void Initial(){
	int i;	
	for(i=0; i<N; i++){
		if(((i-1)<0.5*N)){
			Y[i]=1;
			rho[i]=10.0;
			v[i]=0;
		}else{
			Y[i]=1;
			rho[i]=1.0;
			v[i]=0;
		}
		T[i]=Tdrop+0.5*(1+erf((dx*i-0.5*L)/3))*(Tamp-Tdrop);
		u[i]=rho[i]*Y[i];
		u[i+N]=rho[i]*v[i];
		u[i+N*2]=rho[i]*(CV*T[i]+0.5*v[i]*v[i]);
	}
#if N_SPECIES>1
	int j;
	for(j=1; j<N_SPECIES; j++){ 	
		for(i=0; i<N; i++){
			if(((i-1)<0.5*N)){
				rho[i+j*N]=10.0;
				v[i+j*N]=0;
				Y[i]=1;
			}else{
				rho[i+j*N]=1.0;
				v[i+j*N]=0;
				Y[i]=1;
			}
			T[i+j*N]=Tdrop+0.5*(1+erf((dx*i-0.5*L)/3))*(Tamb-Tdrop);
			u[i+j*3*N]=rho[i+j*N]*Y[i+j*N];
			u[i+N+j*3*N]=rho[i+j*N]*v[i+j*N];
			u[i+N*2+j*3*N]=rho[i+j*N]*(CV*T[i+j*N]+0.5*v[i+j*N]*v[i+j*N]);
		}
	}
#endif
	D[0]=0;
	W[0]=1.0;

}

void Compute_Flux(){
	int i,j; 
	float a;
	float F1, F2, F3;
	float sum1,sum2, sum3;
	
	for(i =1; i<N-1;i++){
		j=0;
		sum1=0;
		sum2=0;
		sum3=0;
		for(j=0; j<N_SPECIES; j++){
			sum1+=W[j]*x[i+j*N];
			//will add compute D_k function here
		}
		for(j=0; j<N_SPECIES; j++){
			sum2+=D[j]*W[j]*(x[i+j*N]-x[i-1+j*N])/dx/sum1;
			sum3+=D[j]*W[j]*CP*T[i+j*N]*(x[i+j*N]-x[i-1+j*N])/dx/sum1;
		}
		for(j=0; j<N_SPECIES; j++){
			F1 = i*dx*i*dx*(u[i+j*3*N]*v[i+j*N]-rho[i+j*N]*D[j]*W[j]*(x[i+j*N]-x[i-1+j*N])/dx/sum1+rho[i+j*N]*Y[i+j*N]*sum2);
			F2 = i*dx*i*dx*(rho[i+j*N]*v[i+j*N]*v[i+j*N] + R*rho[i+j*N]*T[i+j*N] - 1.333*mu*(v[i+j*N]-v[i-1+j*N])/dx);
			F3 = i*dx*i*dx*(v[i+j*N]*(u[i+N*2+j*3*N] + R*rho[i+j*N]*T[i+j*N])-lamda*(T[i+j*N]-T[i-1+j*N])/dx+rho[i+j*N]*sum3);
			fp[i+j*3*N] = 0.5*F1;
			fp[i+N+j*3*N] = 0.5*F2;
			fp[i+N*2+j*3*N] = 0.5*F3;
			fm[i+j*3*N] = 0.5*F1;
			fm[i+N+j*3*N] = 0.5*F2;
			fm[i+N*2+j*3*N] = 0.5*F3;
		}
	}
}

void Compute_U(){
	float FL1, FL2, FL3, FR1, FR2, FR3;
	int i,j,k;
	for(i=1; i<N-1; i++){
		for(j=0; j<N_SPECIES; j++){
			FL1 = fp[i-1+j*3*N]+fm[i+j*3*N];
			FR1 = fp[i+j*3*N]+fm[i+1+j*3*N];
			FL2 = fp[i-1+N+j*3*N]+fm[i+N+j*3*N];
			FR2 = fp[i+N+j*3*N]+fm[i+1+N+j*3*N];
			FL3 = fp[i-1+N*2+j*3*N]+fm[i+N*2+j*3*N];
			FR3 = fp[i+N*2+j*3*N]+fm[i+1+N*2+j*3*N];
			u[i+j*3*N] = u[i+j*3*N] - Z*(FR1-FL1)/i/i/dx/dx;
			u[i+N+j*3*N] = u[i+N+j*3*N] - Z*(FR2-FL2)/i/i/dx/dx;
			u[i+N*2+j*3*N] = u[i+N*2+j*3*N] - Z*(FR3-FL3)/i/i/dx/dx;
		}
	}

	for(j=0; j<N_SPECIES; j++){
		u[0+j*3*N]   = u[1+j*3*N];
		u[N+j*3*N]   = u[1+N+j*3*N];
		u[N*2+j*3*N] = u[1+N*2+j*3*N];
		u[N-1+j*3*N]     = u[N-2+j*3*N];
		u[N-1+N+j*3*N]   = u[N-2+N+j*3*N];
		u[N-1+N*2+j*3*N] = u[N-2+N*2+j*3*N];
	}
	for(i=0; i<N; ++i){
		rho[i+j*N] = 0;
		for(j=0; j<N_SPECIES; j++){
			rho[i+j*N] += u[i+j*3*N];
		}
		for(j=0; j<N_SPECIES; j++){
			Y[i+j*N] = u[i+j*3*N]/rho[i+j*N];
			x[i+j*N] = Y[i+j*N]/W[j];
			v[i+j*N] = u[i+N+j*3*N]/rho[i+j*N];
			T[i+j*N] = (u[i+N*2+j*3*N]/u[i+j*3*N]-v[i+j*N]*v[i+j*N]*0.5)/CV;
		}
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
	free(x);
	free(Y);
	free(T);
	free(W);
}
