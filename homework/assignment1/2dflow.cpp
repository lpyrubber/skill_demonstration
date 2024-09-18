//exercise 5.5
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Minf 0.5
#define Gm 1.4
#define C 1.0
#define TH 0.06
#define Pinf 1.0
#define Rhoinf 1.0 
#define Am (1-Minf*Minf)
#define ainf (sqrt(Gm*Pinf/Rhoinf))
#define Vinf (Minf*ainf)
#define NX 51
#define NY 51
#define NC 21
#define N_it 20
#define N_IT 1000
#define NO (int)(0.5*(NX-NC))
#define L (NX-1)*C
#define D (NY-1)*C
#define LO (L-C)
#define dxmin 0.05*C
#define dymin 0.1*TH
//#define N_METHOD 5
#define N_METHOD 2

void Allocate_Memory();
void Initial();
void Create_Grid();
void Interation();
void Calculate_residual(int method, int time);
void Calculate_residual_old(int method, int time);
void Calculate_CP();
float Calculate_Kappa(int N, float length, float delta);
void Save_Result(int time);
void Free_Memory();

float *x, *y, *phi, *phi_new, *b, *residual;
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
	Calculate_CP();
	fclose(pFile);
	Free_Memory();

}


void Allocate_Memory(){
	x = (float*)malloc(NX * sizeof(float));
	y = (float*)malloc(NY * sizeof(float));
	residual = (float*)malloc(N_METHOD * N_IT * sizeof(float));
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
				phi[i+j*NX+k*NX*NY]=Vinf*x[i];
				phi_new[i+j*NX+k*NX*NY]=Vinf*x[i];
			}
		}
	}
}

void Interation(){
	int i, j, k, l, offset;
	float a, b, c, d;
	for(k=0; k<N_IT; k++){
		//method 1
		offset = 0;
		Calculate_residual(0,k);
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi[i+NX+offset];
			}
		}		
		for(j=1; j<NY-1; j++){
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
				d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
				phi_new[i+j*NX+offset] = 1/(a+b+c+d)*(a*phi[i-1+j*NX+offset]+b*phi[i+1+j*NX+offset]+c*phi[i+(j-1)*NX+offset]+d*phi[i+(j+1)*NX+offset]);
			}
		}
		//method 2
		offset = NX*NY;
		Calculate_residual(1,k);
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi_new[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi_new[i+NX+offset];
			}
		}		
		for(j=1; j<NY-1; j++){
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
				d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
				phi_new[i+j*NX+offset] = 1/(a+b+c+d)*(a*phi_new[i-1+j*NX+offset]+b*phi_new[i+1+j*NX+offset]+c*phi_new[i+(j-1)*NX+offset]+d*phi_new[i+(j+1)*NX+offset]);
			}
		}
		//method 3

		//method 4
		//method 5

		//update phi
		for(l=0; l<N_METHOD; l++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					phi[i+j*NX+l*NX*NY]=phi_new[i+j*NX+l*NX*NY];
				}
			}
		}

	}
}

void Calculate_residual_old(int method, int time){
	int i, j;
	float a, b, c, d, sum, temp, temp1, temp2;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			temp = fabs(phi_new[i+j*NX+method*NX*NY]-phi[i+j*NX+method*NX*NY]);
			sum = (sum>temp) ? sum : temp;
		
		}
	}
	residual[time+method*N_IT] = sum;
}
void Calculate_residual(int method, int time){
	int i, j;
	float a, b, c, d, sum, temp, temp1, temp2;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
			b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
			c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			temp = fabs(-(a+b+c+d)*phi[i+j*NX+method*NX*NY]+a*phi[i-1+j*NX+method*NX*NY]+b*phi[i+1+j*NX+method*NX*NY]+c*phi[i+(j-1)*NX+method*NX*NY]+d*phi[i+(j+1)*NX+method*NX*NY]);
/*
			temp1 = 2*Am*((phi[i+1+j*NX+method*NX*NY]-phi[i+j*NX+method*NX*NY])/(x[i+1]-x[i])\
				     -(phi[i+j*NX+method*NX*NY]-phi[i-1+j*NX+method*NX*NY])/(x[i]-x[i-1]))/(x[i+1]-x[i-1]);
			temp2 = 2*((phi[i+(j+1)*NX+method*NX*NY]-phi[i+j*NX+method*NX*NY])/(y[j+1]-y[j])\
			          -(phi[i+j*NX+method*NX*NY]-phi[i+(j-1)*NX+method*NX*NY])/(y[j]-y[j-1]))/(y[j+1]-y[j-1]);
			temp = fabs(temp1 + temp2);
*/			sum = (sum>temp) ? sum : temp;
		
		}
	}
	residual[time+method*N_IT] = sum;
}

void Calculate_CP(){
	FILE *in2;
	int i,j;
	float u,v,p;
	in2 = fopen("cp.txt","w");
	for(j=0; j<N_METHOD; j++){
		for(i=0;i<NX;i++){
			if(i==0){
				u=(phi[i+1+j*NX*NY]-phi[i+j*NX*NY])/(x[i+1]-x[i]);
			}else if(i==NX-1){
				u=(phi[i+j*NX*NY]-phi[i-1+j*NX*NY])/(x[i]-x[i-1]);

			}else{
				u=(phi[i+1+j*NX*NY]-phi[i-1+j*NX*NY])/(x[i+1]-x[i-1]);
			}
			v=(phi[i+NX+j*NX*NY]-phi[i+j*NX*NY])/(y[1]-y[0]);
			p=Pinf*pow((1-0.5*(Gm-1)*Minf*Minf*((u*u+v*v)/(Vinf*Vinf)-1)),Gm/(Gm-1));
			fprintf(in2,"%e ",2*(Pinf-p)/Rhoinf/Vinf/Vinf);
		}
		fprintf(in2,"\n");
	}
	fclose(in2);

}
void Save_Result(int time){
	int i, j, k;
	FILE* in2;
	if(time == 0){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				fprintf(pFile, "%e ", x[i]);
			}
		}
		fprintf(pFile, "\n");
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				fprintf(pFile, "%e ", y[j]);
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
			}
			fprintf(pFile, "\n");
		}
		in2 = fopen("residual.txt","w");
		for(k=0; k<N_METHOD; k++){
			for(i=0; i<N_IT; i++){
				fprintf(in2, "%e ", residual[i+k*N_METHOD]);
			}
			fprintf(in2,"\n");
		}
		fclose(in2);
	}
}

void Free_Memory(){
	free(phi);
	free(phi_new);
	free(x);
	free(y);
	free(b);
	free(residual);
//	free(A);
}
