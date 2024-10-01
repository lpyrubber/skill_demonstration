#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define LX    100.0
#define R     2.0
#define Theta 1.0

double *x;
double *cx;
double *temp;
int Dim,N_cluster,N_points; 


int main(){
	int i,num,it,j, Nx;
	double temp2;
	FILE *in;
	srand(time(NULL));
	printf("how many dimension you need?\n");
	scanf("%d", &Dim);
	printf("how many cluster you need?\n");
	scanf("%d", &N_cluster);
	printf("how many points you need?\n");
	scanf("%d", &N_points);
	cx = (double*)malloc(N_cluster*Dim*sizeof(double));
	x = (double*)malloc(N_points*Dim*sizeof(double));
	temp = (double*)malloc((Dim+1)*sizeof(double));
	Nx=ceil(N_points/N_cluster);
	for(i=0; i<N_cluster*Dim; i++){
		cx[i]=(double)(LX*rand()/RAND_MAX);
	}
	for(i=0; i<N_cluster; i++){
		it=0;
		num=(i==N_cluster-1)? Nx : N_points-(N_cluster-1)*Nx;
		printf("num=%d,i=%d\n",num, i);
		while(it<num){
			temp[Dim]=0;
			for(j=0;j<Dim;j++){
				temp[j]=double(LX*rand()/RAND_MAX);
				temp[Dim]+=(temp[j]-cx[i+j*N_cluster])*(temp[j]-cx[i+j*N_cluster])/Theta/Theta;
			}
			temp2 = exp(-temp[Dim])/pow(sqrt(M_PI)*Theta,Dim);
			if(rand()/RAND_MAX<temp2){
				for(j=0; j<Dim;j++){
					x[it+i*Nx+j*N_points]=temp[j];
				}
				it++;
			}
		}
	}
	in = fopen("data.txt","w");
	fprintf(in, "%d %d\n",N_points,Dim);
	printf("%d %d\n",N_points,Dim);
	for(i=0; i<N_points; i++){
		printf("%d ",i);
		for(j=0; j<Dim-1; j++){
			fprintf(in, "%lf ",x[i+j*N_points]);
			printf("%lf ",x[i+j*N_points]);
		}
		fprintf(in, "%lf\n",x[i+(Dim-1)*N_points]);
		printf("%lf\n",x[i+(Dim-1)*N_points]);
	}
	fclose(in);
		
	free(cx);
	free(x);
	free(temp);
}


