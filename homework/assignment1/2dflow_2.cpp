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
#define L ((NX-1)*C)
#define D ((NY-1)*C)
#define LO (L-C)
#define dxmin (0.05*C)
#define dymin (0.1*TH)
#define N_METHOD 5

void Allocate_Memory();
void Initial();
void Create_Grid();
void Interation();
void Calculate_residual(int method, int time);
void Calculate_CP();
void LU_Solver(int N);
double Calculate_Kappa(int N, double length, double delta);
void Print_Residual();
void Save_Result(int time);
void Free_Memory();

double *am, *bm, *cm, *dm, *x, *y, *phi, *phi_new, *residual, *cp, *f, *dp, *bp ;
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
	Print_Residual();
	Calculate_CP();
	fclose(pFile);
	Free_Memory();

}


void Allocate_Memory(){
	int N;
	N = (int)((NX+NY+abs(NY-NX))/2);
	printf("N=%d\n",N);
	x = (double*)malloc(NX * sizeof(double));
	y = (double*)malloc(NY * sizeof(double));
	residual = (double*)malloc(N_METHOD * N_IT * sizeof(double));
	phi = (double*)malloc(NX * NY * N_METHOD * sizeof(double));
	phi_new = (double*)malloc(NX * NY * N_METHOD * sizeof(double));
	bp = (double*)malloc(N * sizeof(double));
	cp = (double*)malloc(N * sizeof(double));
	dp = (double*)malloc(N * sizeof(double));
	f = (double*)malloc(N * sizeof(double));
	am = (double*)malloc(N * sizeof(double));
	bm = (double*)malloc(N * sizeof(double));
	cm = (double*)malloc(N * sizeof(double));
	dm = (double*)malloc(N * sizeof(double));
}

void Create_Grid(){
	int i;
	double K;
	for( i = 0; i < NC; i++){
		x[NO + i] = i * dxmin;
	}
	K = Calculate_Kappa(NO+1, LO, dxmin);
	for(i = 0; i < NO; i++){
		x[NO-1-i] = -LO*(exp(K*(i+1)/NO)-1)/(exp(K)-1);
		x[NO + NC + i] = C + LO*(exp(K*(i+1)/NO)-1)/(exp(K)-1); 
	}
	K= Calculate_Kappa(NY, D, dymin);
	y[0]=0;
	for(i = 1; i < NY; i++){
		y[i]=D*(exp(K*(i)/(NY-1))-1)/(exp(K)-1);
	}

}

double Calculate_Kappa(int N, double length, double delta){
	double K=1, i;
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
	double a, b, c, d;
	int flag=1;
	for(k=0; k<N_IT; k++){
		//method 1
		offset = 0;
		Calculate_residual(0,k);
		//left right
		for(j=0; j<NY; j++){
			phi_new[j*NX+offset]=Vinf*x[0];
			phi_new[(NX-1)+j*NX+offset]=Vinf*x[NX-1];
		}
		//top
		for(i=0; i<NX; i++){
			phi_new[i+(NY-1)*NX+offset]=Vinf*x[i];
		}

		for(j=1; j<NY-1; j++){
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
				d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
				phi_new[i+j*NX+offset] = (a*phi[i+1+j*NX+offset]+b*phi[i-1+j*NX+offset]+c*phi[i+(j+1)*NX+offset]+d*phi[i+(j-1)*NX+offset])/(a+b+c+d);
			}
		}
		//bottom
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi[i+NX+offset];
			}
		}

		//method 2
		offset = NX*NY;
		Calculate_residual(1,k);
		//left right
		for(j=0; j<NY; j++){
			phi_new[j*NX+offset]=Vinf*x[0];
			phi_new[(NX-1)+j*NX+offset]=Vinf*x[NX-1];
		}
		//top
		for(i=0; i<NX; i++){
			phi_new[i+(NY-1)*NX+offset]=Vinf*x[i];
		}

		for(j=1; j<NY-1; j++){
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
				d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
				phi_new[i+j*NX+offset] = (a*phi[i+1+j*NX+offset]+b*phi_new[i-1+j*NX+offset]+c*phi[i+(j+1)*NX+offset]+d*phi_new[i+(j-1)*NX+offset])/(a+b+c+d);
			}
		}
		//bottom
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi_new[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi_new[i+NX+offset];
			}
		}

		//method 3
		offset=2*NX*NY;
		Calculate_residual(2,k);
		for(i=1; i<NX-1; i++){
			for(j=0; j<NY; j++){
				//set up tridiagonal matrix
				if(j==0){
					if((i>(NO-1))&&(i<(NO+NC))){
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
					}else{
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=0;
					}
				}else if(j==NY-1){
					am[j]=1;
					bm[j]=0;
					cm[j]=0;
					f[j]=Vinf*x[i];
				}else{
					a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
					c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
					d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
					am[j]=a+b+c+d;
					bm[j]=-c;
					cm[j]=-d;
					f[j]=a*phi[i+1+j*NX+offset]+b*phi[i-1+j*NX+offset];
				}
			}
			LU_Solver(NY);
			for(j=0;j<NY;j++){
				phi_new[i+j*NX+offset]=dm[j];
			}
		}
		//method 4
		offset=3*NX*NY;
		Calculate_residual(3,k);
		for(i=1; i<NX-1; i++){
			for(j=0; j<NY; j++){
				//set up tridiagonal matrix
				if(j==0){
					if((i>(NO-1))&&(i<(NO+NC))){
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
					}else{
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=0;
					}
				}else if(j==NY-1){
					am[j]=1;
					bm[j]=0;
					cm[j]=0;
					f[j]=Vinf*x[i];
				}else{
					a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
					c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
					d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
					am[j]=a+b+c+d;
					bm[j]=-c;
					cm[j]=-d;
					f[j]=a*phi[i+1+j*NX+offset]+b*phi_new[i-1+j*NX+offset];
				}
			}
			LU_Solver(NY);
			for(j=0;j<NY;j++){
				phi_new[i+j*NX+offset]=dm[j];
			}
		}

		//method 5
		offset=4*NX*NY;
		Calculate_residual(4,k);
		for(i=1; i<NX-1; i++){
			for(j=0; j<NY; j++){
				//set up tridiagonal matrix
				if(j==0){
					if((i>(NO-1))&&(i<(NO+NC))){
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
					}else{
						am[j]=1;
						bm[j]=-1;
						cm[j]=0;
						f[j]=0;
					}
				}else if(j==NY-1){
					am[j]=1;
					bm[j]=0;
					cm[j]=0;
					f[j]=Vinf*x[i];
				}else{
					a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
					c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
					d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
					am[j]=a+b+c+d;
					bm[j]=-c;
					cm[j]=-d;
					f[j]=a*phi[i+1+j*NX+offset]+b*phi[i-1+j*NX+offset];
				}
			}
			LU_Solver(NY);
			for(j=0;j<NY;j++){
				phi_new[i+j*NX+offset]=dm[j];
			}
		}

		for(i=0; i<NX; i++){
			phi[i+(NY-1)*NX+offset]=Vinf*x[i];
		}
		for(j=1; j<NY-1; j++){
			for(i=0; i<NX; i++){
				if(i==0||i==NX-1){
					am[i]=1;
					bm[i]=0;
					cm[i]=0;
					f[i]=Vinf*x[i];
				}else{
					a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
					c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
					d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
					am[i]=a+b+c+d;
					bm[i]=-a;
					cm[i]=-b;
					f[i]=c*phi_new[i+(j+1)*NX+offset]+d*phi_new[i+(j-1)*NX+offset];
				}
			}
			LU_Solver(NX);
			for(i=0;i<NX;i++){
				phi[i+j*NX+offset]=dm[i];
			}
		}
		for(i=0; i<NX; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi[i+offset]=phi_new[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi[i+offset]=phi_new[i+NX+offset];
			}
		}
		//update phi
		for(l=0; l<N_METHOD-1; l++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					phi[i+j*NX+l*NX*NY]=phi_new[i+j*NX+l*NX*NY];
				}
			}
		}
		

	}
}

void Calculate_residual(int method, int time){
	int i, j;
	double a, b, c, d, sum, temp, temp1, temp2, temp3,temp4;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			temp1=(phi[i+1+j*NX+method*NX*NY]-phi[i+j*NX+method*NX*NY])/(x[i+1]-x[i]);
			temp2=(phi[i+j*NX+method*NX*NY]-phi[i-1+j*NX+method*NX*NY])/(x[i]-x[i-1]);
			temp3=2*Am*(temp1-temp2)/(x[i+1]-x[i-1]);
			temp1=(phi[i+(j+1)*NX+method*NX*NY]-phi[i+j*NX+method*NX*NY])/(y[j+1]-y[j]);
			temp2=(phi[i+j*NX+method*NX*NY]-phi[i+(j-1)*NX+method*NX*NY])/(y[j]-y[j-1]);
			temp4=2*(temp1-temp2)/(y[j+1]-y[j-1]);
//			sum = (sum > fabs(temp3+temp4)? sum : fabs(temp3+temp4));
			temp = fabs(temp3+temp4);			
/*
			a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
			b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
			c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			temp = fabs(-(a+b+c+d)*phi[i+j*NX+method*NX*NY]+a*phi[i+1+j*NX+method*NX*NY]+b*phi[i-1+j*NX+method*NX*NY]+c*phi[i+(j+1)*NX+method*NX*NY]+d*phi[i+(j-1)*NX+method*NX*NY]);

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
	double u,v,p;
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
				fprintf(in2, "%e ", residual[i+k*N_IT]);
			}
			fprintf(in2,"\n");
		}
		fclose(in2);
	}
}

void LU_Solver(int N){
	int i;
	bp[0]=am[0];
	for(i=1; i<N; i++){
		cp[i]=cm[i]/bp[i-1];
		bp[i]=am[i]-cp[i]*bm[i-1];
	}
	dp[0]=f[0];
	for(i=1; i<N; i++){
		dp[i]=f[i]-cp[i]*dp[i-1];
	}
	dm[N-1]=dp[N-1]/bp[N-1];
	for(i=N-2; i>=0; i--){
		dm[i] = (dp[i]-bm[i]*dm[i+1])/bp[i];
	}

}

void Print_Residual(){
	int i,j,l, cases=6;
	double temp,a,b,c,d;
	int mx[300],my[300],ml[300],index=0;
	FILE *in;

	in=fopen("residual_map.txt", "w");
	for(l=0; l<N_METHOD; l++){
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				if(isnan(phi[i+NX*j+l*NX*NY])){
					printf("phi nan at (x,y) = (%d, %d) with method %d, total index=%d\n",i, j, l, i+NX*j+l*NX*NY);	
				}
				if(i==0||i==NX-1||j==NY-1){
					cases=0;
					temp=fabs(Vinf*x[i]-phi[i+NX*j+l*NX*NY]);
				}else if(j==0){
					if(i>(NO-1) && i<NC+NO){
						cases=1;
						temp=fabs((phi[i+NX+l*NX*NY]-phi[i+l*NX*NY])-Vinf*dymin*((0.5*C-x[i])/sqrt(pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
					}else{
						cases=2;
						temp=fabs(phi[i+NX+l*NX*NY]-phi[i+l*NX*NY]);
					}
				}else{
					cases=3;
                        		a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
                        		c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
                        		d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
                        		temp = fabs(-(a+b+c+d)*phi[i+j*NX+l*NX*NY]+a*phi[i+1+j*NX+l*NX*NY]+b*phi[i-1+j*NX+l*NX*NY]+c*phi[i+(j+1)*NX+l*NX*NY]+d*phi[i+(j-1)*NX+l*NX*NY]);
				}
/*				if(isnan(temp)){
					printf("phi=%e temp nan at (x,y) = (%d, %d), case %d with method %d, total index=%d\n", phi[i+j*NX+l*NX*NY], i, j, cases, l, i+NX*j+l*NX*NY);	
				}
				if(temp>100){
					mx[index]=i;
					my[index]=j;
					ml[index]=l;
					index++;
					printf("%d, %d, %d, %e\n",i,j,l,temp);
				}
*/				fprintf(in, "%e ", temp);
			}
		}
		fprintf(in, "\n");
	}
	fclose(in);
}

void Free_Memory(){
	free(phi);
	free(phi_new);
	free(x);
	free(y);
	free(residual);
	free(cp);
	free(f);
	free(dp);
	free(bp);
	free(am);
	free(bm);
	free(cm);
	free(dm);
}
