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
#define N_METHOD 3

void Allocate_Memory();
void Initial();
void Create_Grid();
void Interation();
void Calculate_residual(int method, int time);
void Calculate_residual_old(int method, int time);
void Calculate_CP();
double Calculate_Kappa(int N, double length, double delta);
void Print_Residual();
void Save_Result(int time);
void Free_Memory();

double *x, *y, *phi, *phi_new, *residual, *cp, *f, *dp ;
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
	cp = (double*)malloc(N * sizeof(double));
	dp = (double*)malloc(N * sizeof(double));
	f = (double*)malloc(N * sizeof(double));
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
//				if((i==0)||(i==NX-1)||(j==0)||(j==NY-1)){
					phi[i+j*NX+k* 
					phi_new[i+j*NX+k* 
//				}else{
//					phi[i+j*NX+k* 
//					phi_new[i+j*NX+k* 
//				}
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
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi[i+NX+offset];
			}
		}		
		for(j=1; j<NY-1; j++){
			c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				phi_new[i+j*NX+offset] = 1/(a+b+c+d)*(a*phi[i-1+j*NX+offset]+b*phi[i+1+j*NX+offset]+c*phi[i+(j-1)*NX+offset]+d*phi[i+(j+1)*NX+offset]);
			}
		}
		//method 2
		offset =  
		Calculate_residual(1,k);
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
				phi_new[i+offset]=phi_new[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				phi_new[i+offset]=phi_new[i+NX+offset];
			}
		}		
		for(j=1; j<NY-1; j++){
			c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				phi_new[i+j*NX+offset] = 1/(a+b+c+d)*(a*phi_new[i-1+j*NX+offset]+b*phi_new[i+1+j*NX+offset]+c*phi_new[i+(j-1)*NX+offset]+d*phi_new[i+(j+1)*NX+offset]);
			}
		}

		//method 3
		offset=2* 
		Calculate_residual(2,k);
		//ydir
		for(i=1; i<NX-1; i++){
			//set up f,
			a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
			b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
			if((i>(NO-1))&&(i<(NC+NO))){
				f[0] = -Vinf*dymin*((0.5*C-x[i])/sqrt((pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));		
			}else{
				f[0]=0;
			}
			for(j=1; j<NY-1; j++){	
				f[j]=-a*phi[i-1+j*NX+offset]-b*phi[i+1+j*NX+offset];
			}
			d=2/(y[NY-1]-y[NY-2])/(y[NY-1]-y[NY-3]);
			f[NY-2] -= d*phi_new[i+(NY-1)*NX+offset];
			//a=c, b=-a-b-c-d, c=d, d=f, 
			cp[0]=-1;
			dp[0]=f[0];
			for(j=1; j<NY-1; j++){
				c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
				d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
				cp[j] = d/(-(a+b+c+d)-c*cp[j-1]);
				dp[j] = (f[j]-c*dp[j-1])/(-(a+b+c+d)-c*cp[j-1]);
			}
			phi_new[i+(NY-2)*NX+offset]=dp[NY-2];
			for(j=NY-3;j>=0;j--){
				phi_new[i+j*NX+offset]=dp[j]-cp[j]*phi_new[i+(j+1)*NX+offset];
			}
		}

/*		//for boundary, ie, j=0;
		for(i=1; i<NX-1; i++){
			if((i>(NO-1))&&(i<(NC+NO))){
		        	phi_new[i+offset]=phi[i+NX+offset]-Vinf*dymin*((0.5*C-x[i])/sqrt(pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2)));
                	}else{
                                phi_new[i+offset]=phi[i+NX+offset];
                        }
		}

		for(j=1; j<NY-1; j++){
			//set up f
			c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
                        d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			
			for(i=1; i<NX-1; i++){
				f[i]=-c*phi[i+(j-1)*NX+offset]-d*phi[i+(j+1)*NX+offset];
				if(isnan(f[i])&&flag){
a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
                        b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
                        c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
                        d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
                        temp = fabs(-(a+b+c+d)*phi[i+j*NX+method* j*NX+method* ]+b i+1 +me  ]+c*phi -1)*NX+ d* ]+d* +( NX+method*  			printf(  d, k=%d\n" k); 
					printf("f=%e\n",f[i]);
					flag=0;
				}	
			}
			a=2*Am/(x[2]-x[0])/(x[1]-x[0]);
                        b=2*Am/(x[NX-1]-x[NX-3])/(x[NX-1]-x[NX-2]);
			f[1]-=a*phi_new[j*NX+offset];
			f[NX-2]-=b*phi_new[NX-1+j*NX+offset];
			
			//a=a, b=-a-b-c-d, c=b, d=f, 
                        b=2*Am/(x[2]-x[0])/(x[2]-x[1]);
			cp[1]=b/(-(a+b+c+d));
			dp[1]=f[1]/(-(a+b+c+d));
			for(i=1; i<NX-1; i++){
				a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				cp[i] = b/(-(a+b+c+d)-b*cp[i-1]);
				dp[i] = (f[i]-a*dp[i-1])/(-(a+b+c+d)-b*cp[i-1]);
				if(isnan(cp[i])||isnan(dp[i])){
					if(flag){
						printf("i=%d, j=%d, k=%d\n",i,j,k);
						printf("cp %e dp %e\n", cp[i],dp[i]);
						printf("%e %e %e %e\n",a,b,c,d);
						printf("%e\n",a+b+c+d+b*cp[i-1]);
						printf("%e\n",f[i]);
						flag=0;
					}
				}
			}
			phi_new[NX-2+j*NX+offset]=dp[NX-2];
			for(i=NX-3;i>0;i--){
				phi_new[i+j*NX+offset]=dp[i]-cp[i]*phi_new[i+1+j*NX+offset];
				if(isnan(phi_new[i])){
					if(flag){
						printf("i=%d, j=%d, k=%d\n",i,j,k);
						printf("cp %e dp %e phi %e\n", cp[i],dp[i],phi_new[i+j*NX+offset]);
						printf("%e %e %e %e\n",a,b,c,d);
						printf("%e\n",a+b+c+d+b*cp[i-1]);
						printf("%e\n",f[i]);
						flag=0;
					}
				}
			}
		}
*/		//method 4
//		offset=3* 
		//method 5
//		offset=4* 
//		if((k&0x1)==0){
			//ydir
//		}else{
			//xdir
		
//		}
		
		
		//update phi
		for(l=0; l<N_METHOD; l++){
			for(j=0; j<NY; j++){
				for(i=0; i<NX; i++){
					phi[i+j*NX+l* j*NX+l* ];   
				}
			}
		}

	}
}

void Calculate_residual_old(int method, int time){
	int i, j;
	double a, b, c, d, sum, temp, temp1, temp2;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			temp = fabs(phi_new[i+j*NX+method* +method* ]);   
			sum = (sum>temp) ? sum : temp;
		
		}
	}
	residual[time+method*N_IT] = sum;
}
void Calculate_residual(int method, int time){
	int i, j;
	double a, b, c, d, sum, temp, temp1, temp2;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
			b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
			c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			temp = fabs(-(a+b+c+d)*phi[i+j*NX+method* j*NX+method* ]+b i+1 +me  ]+c*phi -1)*NX+ d* ]+d* +( NX+method*     
/*
			temp1 = 2*Am*((phi[i+1+j*NX+method* +method* ])/ 1]- \ 
				     -(phi[i+j*NX+method* NX+method* ])/ -x[ )/( ]-x[i-1]);
			temp2 = 2*((phi[i+(j+1)*NX+method* +method* ])/ 1]- \ 
			          -(phi[i+j*NX+method* )*NX+method* ])/ -y[ )/( ]-y[j-1]);
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
				u=(phi[i+1+j* )/ 1]- ; 
			}else if(i==NX-1){
				u=(phi[i+j*  ])/ -x[ ; 

			}else{
				u=(phi[i+1+j*  ])/ 1]- ]); 
			}
			v=(phi[i+NX+j* )/ -y[  
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

void Print_Residual(){
	int i,j,l, cases=6;
	double temp,a,b,c,d;
	int mx[300],my[300],ml[300],index=0;
	FILE *in;

	in=fopen("residual_map.txt", "w");
	for(l=0; l<2; l++){
		for(j=0; j<NY; j++){
                        c=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
                        d=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			for(i=0; i<NX; i++){
				if(isnan(phi[i+NX*j+l* 
					printf("phi nan at (x,y) = (%d, %d) with method %d, total index=%d\n",i, j, l, i+NX*j+l* 
				}
				if(i==0||i==NX-1||j==NY-1){
					cases=0;
					temp=fabs(Vinf*x[i]-phi[i+NX*j+l* 
				}else if(j==0){
					if(i>(NO-1) && i<NC+NO){
						cases=1;
						temp=fabs((phi[i+NX+l* )- dym 0.5 i])/sqrt(pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
					}else{
						cases=2;
						temp=fabs(phi[i+NX+l* );   
					}
				}else{
					cases=3;
					a=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
                        		b=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
                        		temp = fabs(-(a+b+c+d)*phi[i+j*NX+l* j*NX+l* ]+b i+1 +l* phi -1)*NX+ d*phi[i +l* ]);     
				}
				if(isnan(temp)){
					printf("phi=%e temp nan at (x,y) = (%d, %d), case %d with method %d, total index=%d\n", phi[i+j*NX+l* es, l, i+NX*j+l* );	   
				}
				if(temp>100){
					mx[index]=i;
					my[index]=j;
					ml[index]=l;
					index++;
					printf("%d, %d, %d, %e\n",i,j,l,temp);
				}
				fprintf(in, "%e ", temp);
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
}
