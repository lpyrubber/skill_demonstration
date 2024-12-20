#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <omp.h>
#include <math.h>
/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define USE_MATRIX 1
#define N_IT 20
#define SUM_MAX 1e5

void Create_Memory();
char Load_File(char *str);
void Free_Memory();
int Find_Medroid(int tid);
void Label_Point(int tid);
void Save_Result();
static void print_time(double const seconds);
#if USE_MATRIX
void Find_Distance();
#endif

/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/

inline double Calculate_Distance(double *x, int i, int j, int Dim, int N_points){
	double temp=0;
	int k;
	for(k=0; k<Dim; k++){
		temp+=(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
	}
	return sqrtf(temp);
}	

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

int N_c, N_thread, N_points, Dim, flag;
int *label, *c_id, *local_id;
double *x, *sum_dis, *min_c, *local_min;
#if USE_MATRIX
double **distance_m;
#endif
int main(int argc, char** argv){
	int i,j,num;
	double st,et, it, ft;
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
	st=monotonic_seconds();
	it=omp_get_wtime();
//	omp_set_dynamic(0);
	omp_set_num_threads(N_thread);

	for(i=0; i<N_c; i++){
		c_id[i]=i;
		min_c[i]=SUM_MAX;		
	}
	num=N_points/N_thread;
	flag=1;
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		int i=0;
		int temp=0;
#if USE_MATRIX
		Find_Distance();
#endif
//		int N_local=(tid<N_points-num*N_thread)?num+1:num;
//		printf("tid=%d, N_local=%d\n",tid, N_local);

//		local_min=(double*)malloc(N_c*sizeof(double));
//		local_id=(int*)malloc(N_c*sizeof(double));
		Label_Point(tid);
		while(flag && (i<N_IT)){
			temp=Find_Medroid(tid);
//			printf("temp=%d at %d\n",temp, tid);
			#pragma omp barrier
			#pragma omp critical
			{
				flag&=temp;
//				printf("flag=%d at %d\n",flag, tid);
			}
			#pragma omp barrier
			#pragma omp master
			{
				flag=!flag;
			}
			#pragma omp barrier
//			printf("final flag=%d at %d\n",flag, tid);
			Label_Point(tid);
			i++;
		}
		if(tid==0)	printf("break at it=%d\n",i-1);
//		free(local_id);
//		free(local_min);

	}
	et=monotonic_seconds();
	ft=omp_get_wtime();
	printf("\n\nTime taken is %f\n",ft-it);
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
			if(fscanf(in, "%f ", x+(i+j*N_points))<1){
				printf("can't get data\n");
				return 3;
			}
		}
		if(fscanf(in, "%f\n", x+(i+(Dim-1)*N_points))<1){
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
	local_min=(double*)malloc(N_thread*N_c*sizeof(double));
	label=(int*)malloc(N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
	local_id=(int*)malloc(N_thread*N_c*sizeof(int));
#if USE_MATRIX
	int i;
	distance_m=(double**)malloc(N_points*sizeof(double*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(double*)malloc(N_points*sizeof(double));
	}
#endif

}
#if USE_MATRIX
void Find_Distance(){
	int i, j, k;
	#pragma omp for //schedule(dynamic)
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
void Label_Point(int tid){
	int i, j, k;
	double sum;
	double temp;
	#pragma omp barrier
	#pragma omp for // schedule(dynamic)
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

int Find_Medroid(int tid){
	int i,j,k;
	int temp;
	int flag;
	
	flag=1;
	for(i=0; i<N_c; i++){
		local_min[i+tid*N_c]=min_c[i];
		local_id[i+tid*N_c]=c_id[i];
//		printf("at %d,cluster %i,min=%f,id=%d\n",tid,i,local_min[i+tid*N_c],local_id[i+tid*N_c]);
	}
	#pragma omp for // schedule(dynamic)
	for(i=0; i<N_points; i++){
		sum_dis[i]=0;
	}
	#pragma omp barrier
	#pragma omp for// schedule(dynamic)
	for(i=0; i<N_points; i++){
		for(j=0; j<N_points; j++){
			if(label[i]==label[j]){
#if USE_MATRIX
				sum_dis[i]+=distance_m[i][j]/N_points*N_c;
#else
				temp=Calculate_Distance(x,i,j,Dim,N_points);		
				sum_dis[i]+=temp;
#endif				
			}
		}
		if(sum_dis[i]<local_min[label[i]+tid*N_c]){
			local_id[label[i]+tid*N_c]=i;
			local_min[label[i]+tid*N_c]=sum_dis[i];
		}
	}
	#pragma omp barrier
	#pragma omp critical
	{
		for(i=0; i<N_c; i++){
			if(local_min[i+tid*N_c]<min_c[i]){
				c_id[i]=local_id[i+tid*N_c];
				min_c[i]=local_min[i+tid*N_c];
				flag=0;
//				printf("flag=%d\n at tid=%d\n",flag,tid);
			}
		}

	}
		
	return flag;
}

void Save_Result(){
	FILE *out;
	int i, j;
	out = fopen("label.txt","w");
	for(i=0; i<N_points; i++){
		fprintf(out, "%d\n",label[i]);
	}
	fclose(out);
	out = fopen("medroid.txt","w");
	for(i=0; i<N_c; i++){
		for(j=0; j<Dim-1;j++){
			fprintf(out, "%f ", x[c_id[i]+j*N_points]);
		}
		fprintf(out, "%f\n",x[c_id[i]+(Dim-1)*N_points]);
	}
	fclose(out);
}


void Free_Memory(){
	free(x);
	free(sum_dis);
	free(min_c);
	free(label);
	free(c_id);
	free(local_id);
	free(local_min);
#if USE_MATRIX
	int i;
	for(i=0; i<N_points; i++){
		free(distance_m[i]);
	}
	free(distance_m);
#endif
}
