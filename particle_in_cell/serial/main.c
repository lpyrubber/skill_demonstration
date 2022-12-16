#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#define T 100 
#define dt 0.01
#define Nt ((int)(T/dt))
#define N 1e4
#define L M_PI*4
#define Nx 200
#define dx ((double)(L/Nx))
#define nppc ((double)(N/Nx))
#define kt 10 
#define Itstep 1e9
#define Itlimit 1e-6

void Memset();
void Initial();
void Compute();
void Free();
void Pushx();
void Pushv();
void Hist();
void Efield();
void Save();

double *x, *v, *a, *e, *rho;

int main(){
	double a;
	a=(double)N;
	printf("%f\n",a);
	a=(double)Nt;
	printf("%f\n",a);
	Memset();
	Initial();
	Compute();
	Free();
	return 0;
}

void Initial(){
	srand(time(NULL));
	double a,t,f;
	int i,j;
	for(i=0; i<N; i++){
		//for landau
/*		for(j=1;j>0;j++){
			a=5-10*(double)rand()/RAND_MAX;
			f=exp(-(a)*(a)*0.5)/sqrt(M_PI/2);			
			t=(double)rand()/RAND_MAX;
			if(f>=t){
				//printf("%e\n",t);
				v[i]=a;
				break;	
			}
		}
*/		//for two stream
		if(i%2==0){
				v[i]=1;
		}else{
				v[i]=-1;
		}
		for(j=1;j>0;j++){
			a=L*(double)rand()/RAND_MAX;
			f=1+0.5*cos(0.5*a);
			t=(double)rand()/RAND_MAX;
			if(f>=t){
				x[i]=a;
				break;	
			}
		}
		
	}
	Save();
	printf("initializaton complete\n");
}

void Memset(){
	x=(double*)malloc(N*sizeof(double));
	v=(double*)malloc(N*sizeof(double));
	a=(double*)malloc(N*sizeof(double));
	rho=(double*)malloc(Nx*sizeof(double));
	e=(double*)malloc(Nx*sizeof(double));
	printf("Allocation complete\n");
}

void Compute(){
	int i;
	for(i=0; i<Nt; i++){
		Pushx();
		Hist();
		Efield();
		Pushv();
		Pushx();
		if(i%kt==(kt-1)){
			Save();	
		}
		if(i%50==0){
			printf("finish %2.2f percent \n",(double)(i)/Nt*100.0);
		}
	}
	printf("computation complete\n");
}

void Pushx(){
	int i;
	for(i=0;i<N;i++){
		x[i] = x[i]+v[i]*dt*0.5;
		if(x[i]<0){
			x[i]=x[i]+L;
		}
		if(x[i]>L){
			x[i]=x[i]-L;
		}
	}
}

void Pushv(){
	int i,j;
	for(i=0;i<N;i++){
		j=(int)(x[i]/dx);
		v[i]=v[i]-e[j]*dt;
	}
}

void Efield(){	
	//CG method
	double *r, *p, *q, *phi, b;
	double sum, sum1, sum2, sum3, alfa, beta, temp, error;
	int k,i, up, down;
	r = (double*)malloc(sizeof(double)*Nx);
	p = (double*)malloc(sizeof(double)*Nx);
	q = (double*)malloc(sizeof(double)*Nx);
	phi = (double*)malloc(sizeof(double)*Nx);
	for(i=0; i<Nx; i++){
		phi[i]=0;
	}
	for(i=0; i<Nx; i++){
		up=(i+Nx+1)%Nx;
		down=(i+Nx-1)%Nx;
		b=(rho[i]-1)*dx*dx;
		r[i]=b-phi[down]-phi[up]+2*phi[i];
		p[i]=r[i];
	}
	for(i=0; i<Nx; i++){
		up=(i+Nx+1)%Nx;
		down=(i+Nx-1)%Nx;
		q[i]=r[up]+r[down]-2*r[i];
	}
	sum1=0;
	sum2-0;
	for(i=0; i<Nx; i++){
		sum1+=r[i]*r[i];
		sum2+=p[i]*q[i];
	}
	alfa=sum1/sum2;
	for(k=0; k<Itstep; k++){
		sum=0;
		for(i=0; i<Nx; i++){
			temp=phi[i];
			phi[i] = phi[i]+alfa*p[i];
			temp-=phi[i];
			sum+=temp*temp;
			r[i] = r[i]-alfa*q[i];
		}
		sum3=0;
		for(i=0; i<Nx; i++){
			sum3+=r[i]*r[i];
		}
		beta=sum3/sum1;
		for(i=0; i<Nx; i++){
			p[i] = r[i] + beta*p[i];
		}
		sum1=sum3;
		sum2=0;
		for(i=0; i<Nx; i++){
			up=(i+Nx+1)%Nx;
			down=(i+Nx-1)%Nx;
			q[i] = p[down]+p[up]-2*p[i];
			sum2+=q[i]*p[i];
		}
		alfa=sum1/sum2;
		error = sqrtf(sum);
		if(error<Itlimit){
			break;
		}
	}
	for(i=0;i<Nx;i++){
		up=(i+Nx+1)%Nx;
		down=(i+Nx-1)%Nx;
		e[i]=-0.5*(phi[up]-phi[down])/dx;
	}
	free(r);
	free(p);
	free(q);
	free(phi);
}

void Hist(){
	int i;
	int up,down;
	for(i=0; i<Nx; i++){
		rho[i]=0;
	}
	for(i=0; i<N; i++){
		if((x[i]<0.5*dx)||(x[i]>(L-0.5*dx))){
			up=0;
			down=Nx-1;
		}else{
			up=(int)(x[i]/dx-0.5)+1;
			down=up-1;
		}
		if(x[i]<0.5*dx){
			rho[up]+=(0.5*dx+x[i])/dx/nppc;
			rho[down]+=(0.5*dx-x[i])/dx/nppc;
		}else if(x[i]>(L-0.5*dx)){
			rho[up]+=(x[i]-L+0.5*dx)/dx/nppc;
			rho[down]+=(0.5*dx+L-x[i])/dx/nppc;
		}else{
			rho[up]+=(x[i]-down*dx-0.5*dx)/dx/nppc;
			rho[down]+=(up*dx+0.5*dx-x[i])/dx/nppc;
		}
	}
}

void Save(){
	static int i=0;	
	FILE *in,*inn;
	int j;
	char buffer[20];
//	printf("save file\n");
	if(i==0){
		sprintf(buffer,"init_data.txt");
	}else{
		sprintf(buffer,"data_%d.txt",i);
	}
	in = fopen(buffer,"w");
	for(j=0;j<N;j++){
		fprintf(in,"%f\t%f\n",x[j],v[j]);
	}
	if(i>=1){
		sprintf(buffer,"e_data_%d.txt",i);
		inn = fopen(buffer,"w");
		for(j=0;j<Nx;j++){
			fprintf(inn,"%f\t%f\n",j*dx+0.5*dx,e[j]);
		}
		fclose(inn);
	}
	fclose(in);
	i+=1;
}

void Free(){
	free(x);
	free(v);
	free(a);
	free(e);
	free(rho);
}
