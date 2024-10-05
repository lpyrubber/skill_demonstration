#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <math.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define N_IT 20
#define SUM_MAX 1e5
#define USE_MATRIX 1

void Create_Memory();
char Load_File(char *str);
void Free_Memory();
int Find_Medroid();
void Label_Point();
void Save_Result();
static void print_time(double const seconds);
#if USE_MATRIX
void Find_Distance();
#endif

inline double Calculate_Distance(double *x, int i, int j, int Dim, int N_points){
	double temp=0;
	int k;
	for(k=0; k<Dim; k++){
		temp+=(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
	}
	return sqrtf(temp);
}

/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/


static inline double monotonic_seconds()
{
#ifdef __MACH__
  /* OSX */
  static mach_timebase_info_data_t info;
  static double seconds_per_unit;
  if(seconds_per_unit == 0) {
    mach_timebase_info(&info);
    seconds_per_unit = (info.numer / info.denom) / 1e9;
  }
  return seconds_per_unit * mach_absolute_time();
#else
  /* Linux systems */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

int N_c, N_thread, N_points, Dim;
int *label, *c_id;
double *x, *sum_dis, *min_c;
#if USE_MATRIX
float **distance_m;
#endif
int main(int argc, char** argv){
	int i,j,k, flag;
	double st,et;
	if(argc<4) {
		printf("Not enough of arguemt\n");
		return 1;
	}
	if(argc>4) {
		printf("Too many arguements\n");
		return 2;
	}
	
	N_c=atoi(argv[2]);
	N_thread = atoi(argv[3]);
	if(Load_File(argv[1])){
		return 3;
	}
	printf("N_c =%d, N_thread = %d\n", N_c, N_thread);
	//intialization
	st=monotonic_seconds();
#if USE_MATRIX
		Find_Distance();
#endif
	for(i=0; i<N_c; i++){
		c_id[i]=i;
		min_c[i]=SUM_MAX;
	}
	Label_Point();
	i=0;
	flag=1;
	while(flag&&i<N_IT){
		flag=Find_Medroid();
		flag=1;
		Label_Point();
		i++;
	}
	printf("break at %d\n",i-1);
	et=monotonic_seconds();
	print_time(et-st);
	Save_Result();
	Free_Memory();
	return 0;
}




/**
* @brief Output the seconds elapsed while clustering.
*
* @param seconds Seconds spent on k-means clustering, excluding IO.
*/
static void print_time(double const seconds)
{
  printf("k-means clustering time: %0.04fs\n", seconds);
}

char Load_File(char *str){
	FILE *in;
	int i,j;
	in = fopen(str, "r");
	if(in == NULL){
		printf("can't find the file\n");
		return 1;
	}
	if(fscanf(in,"%d %d\n",&N_points, &Dim )==2){
		printf("N_point = %d, Dim = %d\n",N_points, Dim);
	}else{
		printf("can't get Number of nodes and relative dimension\n");
		return 2;
	}
	Create_Memory();
	for(i=0; i<N_points; i++){
		for(j=0; j<Dim-1; j++){
			if(fscanf(in, "%lf ", x+(i+j*N_points))<1){
				printf("can't get data\n");
				return 3;
			}
		}
		if(fscanf(in, "%lf\n", x+(i+(Dim-1)*N_points))<1){
			printf("can't get data\n");
			return 3;
		}
	}
		
	return 0;
}


void Create_Memory(){
	x=(double*)malloc(N_points*Dim*sizeof(double));
	sum_dis=(double*)malloc(N_points*sizeof(double));
	min_c=(double*)malloc(N_c*sizeof(double));
	label=(int*)malloc(2*N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
#if USE_MATRIX
	int i;
	distance_m=(float**)malloc(N_points*sizeof(float*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(float*)malloc(N_points*sizeof(float));
	}
#endif	
}
#if USE_MATRIX
void Find_Distance(){
	int i, j, k;
	#pragma omp for
	for(i=0; i<N_points; i++){
		for(j=i; j<N_points; j++){
			distance_m[i][j]=0;
			for(k=0; k<Dim; k++){
				distance_m[i][j]+=(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
			}
			distance_m[i][j]=sqrtf(distance_m[i][j]);
			distance_m[j][i]=distance_m[i][j];
		}
	}
}
#endif
void Label_Point(){
	int i, j, k;
	double sum,temp;
	for(i=0; i<N_points; i++){
		sum=SUM_MAX;
		label[i]=-1;
		for(j=0; j<N_c; j++){
#if USE_MATRIX
			if(distance_m[i][c_id[j]]<sum){
				label[i]=j;
				sum=distance_m[i][c_id[j]];
			}

#else			
			temp=Calculate_Distance(x,i,c_id[j],Dim,N_points);
			if(temp<sum){
				label[i]=j;
				sum=temp;
			}
#endif
		}	
	}

}

int Find_Medroid(){
	int flag=0;
	int i,j,k;
	double temp=0;
	for(i=0; i<N_points; i++){
		sum_dis[i]=0;
		for(j=0; j<N_points; j++){
			if(label[i]==label[j]){
#if USE_MATRIX
				sum_dis[i]+=distance_m[i][j];
#else
				temp=Calculate_Distance(x,i,j,Dim,N_points);		
				sum_dis[i]+=temp;
#endif						

			}
		}
		if(sum_dis[i]<min_c[label[i]]){
			c_id[label[i]]=i;
			min_c[label[i]]=sum_dis[i];
			flag++;
		}
	}
	return (flag>0)?1:0;
}

void Save_Result(){
	FILE *out;
	int i,j;
	out = fopen("label.txt","w");
	for(i=0; i<N_points; i++){
		fprintf(out, "%d\n",label[i]);
	}
	fclose(out);
	out = fopen("medroid.txt","w");
	for(i=0; i<N_c; i++){
		for(j=0; j<Dim-1; j++){
			fprintf(out, "%lf ",x[c_id[i]+j*N_points]);
		}
		fprintf(out,"%lf\n",x[c_id[i]+(Dim-1)*N_points]);
	}
	fclose(out);
}


void Free_Memory(){
	free(x);
	free(sum_dis);
	free(min_c);
	free(label);
	free(c_id);
#if USE_MATRIX
	int i;
	for(i=0; i<N_points; i++){
		free(distance_m[i]);
	}
	free(distance_m);
#endif
}
