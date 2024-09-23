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


void Allocate_Memory();
void Initial();
void Create_Grid();
void Interation();
void Calculate_residual(int time);
void Calculate_CP();
void LU_Solver(int N);
double Calculate_Kappa(int N, double length, double delta);
void Print_Residual();
void Save_Result(int time);
void Free_Memory();

double *am, *bm, *cm, *dm, *x, *y, *phi, *phi_new, *residual, *cp, *f, *dp, *bp;
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
	residual = (double*)malloc(N_IT * sizeof(double));
	phi = (double*)malloc(NX * NY * sizeof(double));
	phi_new = (double*)malloc(NX * NY * sizeof(double));
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
	
	for(j=0; j<NY; j++){
		for(i=0; i<NX; i++){
				phi[i+j*NX]=1.0;
				phi_new[i+j*NX]=1.0;
		}
	}
	
}

void Interation(){
	int i, j, k, l;
	double a1, b1, c1, d1, a2, b2, u, v, Am1d, Am2d, phi1d, phi2d, mu1d, mu2d;
	for(k=0; k<N_IT; k++){
		Calculate_residual(k);
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
					f[j]=1.0;
				}else{
					phi1d=(phi[i+1+j*NX]-phi[i-1+j*NX])/(x[i+1]-x[i-1]);
					Am1d=1-Minf*Minf-phi1d*(Gm+1)*Minf*Minf/Vinf;
					mu1d=(Am1d>0)?0:1;
					a1=2/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
					b1=2/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
					c1=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
					d1=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);

					if(i==1){
						am[j]=-Am1d*(a1+b1)-c1-d1;
						bm[j]=c1;
						cm[j]=d1;
						f[j]=-Am1d*a1*phi[i+1+j*NX]-Am1d*b1*phi[i-1+j*NX];
					}else{
						phi2d=(phi[i+j*NX]-phi[i-2+j*NX])/(x[i]-x[i-2]);
						Am2d=1-Minf*Minf-phi2d*(Gm+1)*Minf*Minf/Vinf;
						mu2d=(Am2d>0)?0:1;
						a2=2/(x[i]-x[i-2])/(x[i]-x[i-1]);
						b2=2/(x[i]-x[i-2])/(x[i-1]-x[i-2]);

						am[j]=-((1-mu1d)*Am1d*(a1+b1)-mu2d*Am2d*a2+c1+d1);
						bm[j]=c1;
						cm[j]=d1;
						f[j]=-((mu2d*Am2d*b2)*phi[i-2+j*NX]+((1-mu1d)*Am1d*b1-mu2d*Am2d*(a2+b2))*phi[i-1+j*NX]+((1-mu1d)*Am1d*a1)*phi[i+1+j*NX]);
					}
				}
			}
			LU_Solver(NY);
			for(j=0;j<NY;j++){
				phi_new[i+j*NX]=dm[j];
			}
		}
		for(j=0; j<NY; j++){
			for(i=0; i<NX; i++){
				phi[i+j*NX]=phi_new[i+j*NX];
			}
		}
		

	}
}

void Calculate_residual(int time){
	int i, j;
	double a1, b1, c1, d1, a2, b2, Am1d, Am2d, mu1d, mu2d, sum, temp;
	double phi1d, phi2d;
	sum=0;
	for(i=1; i<NX-1; i++){
		for(j=1; j<NY-1; j++){
			phi1d=(phi[i+1+j*NX]-phi[i-1+j*NX])/(x[i+1]-x[i-1]);
			Am1d=1-Minf*Minf-phi1d*(Gm+1)*Minf*Minf/Vinf;
			mu1d=(Am1d>0)?0:1;
			a1=2/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
			b1=2/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
			c1=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
			d1=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
			if(i==1){
				temp=(-Am1d*(a1+b1)-c1-d1)*phi[i+j*NX]+c1*phi[i+(j+1)*NX]+d1*phi[i+(j-1)*NX]-(-Am1d*a1*phi[i+1+j*NX]-Am1d*b1*phi[i-1+j*NX]);
			}else{
				phi2d=(phi[i+j*NX]-phi[i-2+j*NX])/(x[i]-x[i-2]);
				Am2d=1-Minf*Minf-phi2d*(Gm+1)*Minf*Minf/Vinf;
				mu2d=(Am2d>0)?0:1;
				a2=2/(x[i]-x[i-2])/(x[i]-x[i-1]);
				b2=2/(x[i]-x[i-2])/(x[i-1]-x[i-2]);
				temp=-((1-mu1d)*Am1d*(a1+b1)-mu2d*Am2d*a2+c1+d1)*phi[i+j*NX]+c1*phi[i+(j+1)*NX]+d1*phi[i+(j-1)*NX]+((mu2d*Am2d*b2)*phi[i-2+j*NX]+((1-mu1d)*Am1d*b1-mu2d*Am2d*(a2+b2))*phi[i-1+j*NX]+((1-mu1d)*Am1d*a1)*phi[i+1+j*NX]);
			}
			temp = fabs(temp);
			sum = (sum > temp) ? sum : temp;
		}
	}
	residual[time] = sum;
}

void Calculate_CP(){
	FILE *in2;
	int i,j;
	double u,v,p;
	in2 = fopen("cp.txt","w");
	for(i=0;i<NX;i++){
		if(i==0){
			u=(phi[i+1]-phi[i])/(x[i+1]-x[i]);
		}else if(i==NX-1){
			u=(phi[i]-phi[i-1])/(x[i]-x[i-1]);
		}else{
			u=(phi[i+1]-phi[i-1])/(x[i+1]-x[i-1]);
		}
		v=(phi[i+NX]-phi[i])/(y[1]-y[0]);
		p=Pinf*pow((1-0.5*(Gm-1)*Minf*Minf*((u*u+v*v)/(Vinf*Vinf)-1)),Gm/(Gm-1));
		fprintf(in2,"%e ",2*(Pinf-p)/Rhoinf/Vinf/Vinf);
	}
	fprintf(in2,"\n");
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
		
		for(j = 0; j < NY; j++){
			for( i = 0; i < NX; i++){
				fprintf(pFile, "%e ", phi[ i + j * NX]);
			}
		}
		fprintf(pFile, "\n");
		in2 = fopen("residual.txt","w");
		
		for(i=0; i<N_IT; i++){
			fprintf(in2, "%e ", residual[i]);
		}
		fprintf(in2,"\n");
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
/*	int i,j,l, cases=6;
	double temp,a,b,c,d;
	FILE *in;

	in=fopen("residual_map.txt", "w");
	for(j=0; j<NY; j++){
		for(i=0; i<NX; i++){
			if(i==0||i==NX-1||j==NY-1){
				cases=0;
				temp=fabs(Vinf*x[i]-phi[i+NX*j]);
			}else if(j==0){
				if(i>(NO-1) && i<NC+NO){
					cases=1;
					temp=fabs((phi[i+NX]-phi[i])-Vinf*dymin*((0.5*C-x[i])/sqrt(pow(0.25*C*C/TH+0.25*TH,2)-pow(x[i]-0.5*C,2))));
				}else{
					cases=2;
					temp=fabs(phi[i+NX]-phi[i]);
				}
			}else{
				cases=3;
                       		a=2*Am/(x[i+1]-x[i-1])/(x[i+1]-x[i]);
				b=2*Am/(x[i+1]-x[i-1])/(x[i]-x[i-1]);
                       		c=2/(y[j+1]-y[j-1])/(y[j+1]-y[j]);
                       		d=2/(y[j+1]-y[j-1])/(y[j]-y[j-1]);
                       		temp = fabs(-(a+b+c+d)*phi[i+j*NX]+a*phi[i+1+j*NX]+b*phi[i-1+j*NX]+c*phi[i+(j+1)*NX]+d*phi[i+(j-1)*NX]);
			}
			fprintf(in, "%e ", temp);
		}
	}
	fprintf(in, "\n");
	fclose(in);
*/
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
