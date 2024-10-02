#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define LX    100.0
#define Theta 0.5

double *x;
double *cx;
double *temp;
int Dim,N_cluster,N_points; 


int main(){
	int i,num,it,j, Nx, offset;
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
	Nx=(int)(N_points/N_cluster);
	printf("Nx=%d\n",Nx);
	for(i=0; i<N_cluster*Dim; i++){
		cx[i]=(double)(LX*rand()/RAND_MAX);
	}
	offset=0;
	for(i=0; i<N_cluster; i++){
		it=0;
		num=(i<(N_points-Nx*N_cluster))?Nx+1 : Nx;
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
					x[it+offset+j*N_points]=temp[j];
				}
				it++;
			}
		}
		offset+=num;
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


