#include <stdio.h>
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
#define NO_STEP  100
#define mu 1.0
#define lamda 1.0

float *u, *rho, *v, *T, *D, *Y, *x, *W;
float *fp;
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
	int i;
	Allocate_Memory();
	Initial();
	Save_Result();
	for(i=0; i<NO_STEP; i++){
		Compute_Flux();
		Compute_U();
		if(i%100==99){
		       Save_Result();
		}	       
	}
	printf("%f %f\n", u[0], u[N]);
	printf("%f %f\n", fp[1], fp[2]);
	Save_Result();
	Free_Memory();
	return 0;
}

void Allocate_Memory(){
	size_t size=N*sizeof(float);

	rho = (float*)malloc(size);
	v = (float*)malloc(size);
	Y = (float*)malloc(N_SPECIES*size);
	x = (float*)malloc(N_SPECIES*size);
	T = (float*)malloc(size);
#if N_SPECIES>1
	u = (float*)malloc((3+N_SPECIES)*size);
	fp = (float*)malloc((3+N_SPECIES)*size);
	printf("u and fp have %d rows\n", 3+N_SPECIES);
#else
	u = (float*)malloc(3*size);
	fp = (float*)malloc(3*size);
	printf("u and fp have %d rows\n", 3);
#endif
	D = (float*)malloc(N_SPECIES*sizeof(float));
	W = (float*)malloc(N_SPECIES*sizeof(float));
}

void Initial(){
	int i,j;
	float W_total;
        D[0]=0;
	W[0]=1.0;
	W_total=0.0;
#if N_SPECIES>1	
	//calculate for W
	for(j=0; j<N_SPECIES; i++){
		W_total+=1/N_SPECIES/W[j];
	}
	printf("calculate W\n");
#endif
	for(i=0; i<N; i++){
		if(((i-1)<0.5*N)){
			rho[i]=10.0;
			v[i]=0;
		}else{
			rho[i]=1.0;
			v[i]=0;
		}
		T[i]=Tdrop+0.5*(1+erf((dx*i-0.5*L)*5))*(Tamp-Tdrop);
		u[i]=rho[i];
		u[i+N]=rho[i]*v[i];
		u[i+2*N]=rho[i]*(CV*T[i]+0.5*v[i]*v[i]);
#if N_SPECIES>1
		printf("calculating Y,x\n");
		for(j=0;j<N_SPECIES;j++){
			Y[i+j*N]=1/N_SPECIES;
			x[i+j*N]=Y[i+j*N]/W[j]/W_total;
			u[i+(j+3)*N]=rho[i]*Y[i+j*N]
		}
#endif
	}
}

void Compute_Flux(){
	int i,j; 
	float F1, F2, F3;
	float sum1,sum2, sum3;
	
	for(i=1; i<N;i++){
#if N_SPECIES>1
		sum1=0;
		sum2=0;
		sum3=0;
		for(j=0; j<N_SPECIES; j++){
			sum1+=W[j]*x[i+j*N];
			//will add compute D_k function here
		}
		for(j=0; j<N_SPECIES; j++){
			sum2+=D[j]*W[j]*(x[i+j*N]-x[i-1+j*N])/dx/sum1;
			sum3+=D[j]*W[j]*CP*T[i]*(x[i+j*N]-x[i-1+j*N])/dx/sum1;
		}
#endif
		F1 = i*dx*i*dx*(rho[i]*v[i]); 
		//pressure term might have mistake since 
		F2 = i*dx*i*dx*(rho[i]*v[i]*v[i] + R*rho[i]*T[i] - 1.333*mu*(v[i]-v[i-1])/dx);
		F3 = i*dx*i*dx*(v[i]*(u[i] + R*rho[i]*T[i])-lamda*(T[i]-T[i-1])/dx);
#if N_SPECIES>1
		F3 -=i*dx*i*dx*(rho[i]*sum3)
#endif
		fp[i] = 0.5*F1;
		fp[i+N] = 0.5*F2;
		fp[i+2*N] = 0.5*F3;
#if N_SPECIES>1
		for(j=0; j<N_SPECIES; j++){
			F1 = i*dx*i*dx*(u[i+(3+j)*N]*v[i]-rho[i]*D[j]*W[j]*(x[i+j*N]-x[i-1+j*N])/dx/sum1+rho[i]*Y[i+j*N]*sum2);
			fp[i+(j+3)*N] = 0.5*F1;
		}
#endif
	}

}

void Compute_U(){
	float FL1, FL2, FL3, FR1, FR2, FR3, W_total;
	int i,j,k;
	for(i=2; i<N-1; i++){
		FL1 = fp[i-1]+fp[i];
		FR1 = fp[i]+fp[i+1];
		FL2 = fp[i-1+N]+fp[i+N];
		FR2 = fp[i+N]+fp[i+1+N];
		FL3 = fp[i-1+N*2]+fp[i+N*2];
		FR3 = fp[i+N*2]+fp[i+1+N*2];
		u[i] = u[i] - Z*(FR1-FL1)/i/i/dx/dx;
		u[i+N] = u[i+N] - Z*(FR2-FL2)/i/i/dx/dx;
		u[i+N*2] = u[i+N*2] - Z*(FR3-FL3)/i/i/dx/dx;
#if N_SPECIES>1
		for(j=0;j<N_SPECIES;j++){
			FL1 = fp[i-1+(j+3)*N]+fp[i+(j+3)*N];
			FR1 = fp[i+(j+3)*N]+fp[i+1+(j+3)*N];
			u[i+(j+3)*N] = u[i+(j+3)*N]+Z*(FR1-FL1)/i/i/dx/dx;
		}
#endif
	}
	u[1]     = u[2];
	u[1+N]   = u[2+N];
	u[1+N*2] = u[2+N*2];	
	u[0]   = u[1];
	u[N]   = u[1+N];
	u[N*2] = u[1+N*2];
	u[N-1]     = u[N-2];
	u[N-1+N]   = u[N-2+N];
	u[N-1+N*2] = u[N-2+N*2];
#if N_SPECIES>1
	for(j=0; j<N_SPECIES; j++){
		u[N*(j+3)] = u[1+N*(j+3)];
		u[N-1+N*(j+3)] = u[N-2+N*(j+3)];
	}
#endif
	for(i=0; i<N; i++){
		rho[i] = u[i];
		v[i] = u[i+N]/rho[i];
		T[i] = (u[i+N*2]/u[i]-v[i]*v[i]*0.5)/CV;
#if N_SPECIES>1
		W_total=0;
		for(j=0; j<N_SPECIES; j++){
			Y[i+j*N] = u[i+(j+3)*N]/rho[i];
			W_total+=Y[i+j*N]/W[j];
		}
		for(j=0; j<N_SPECIES; j++){
			x[i+j*N]=Y[i+j*N]/W[j]/W_total;
		}
#endif
	}
}

void Save_Result(){
	int i,j;
	static int time=0;
  	char temp[20];
        sprintf(temp, "data_%d.txt",time);
        pFile = fopen(temp, "w");
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=0; i<N; i++){
			fprintf(pFile, "%e %e %e %e", i*dx, rho[i], v[i], T[i]);
#if N_SPECIES>1
			for(j=0; j<N_SPECIES; j++){
				fprintf(pFile, " %e %e", Y[i+j*N], x[i+j*N]);
			}
#endif
			fprintf(pFile, "\n");
		}
	}
	fclose(pFile);
	time++;
}

void Free_Memory(){
	free(fp);
	free(u);
	free(rho);
	free(v);
	free(x);
	free(Y);
	free(T);
	free(W);
}
