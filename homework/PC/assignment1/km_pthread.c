#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <pthread.h>
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
int Find_Medroid(int tid, int num, int offset);
void Label_Point(int tid, int num, int offset);
void Save_Result();
static void print_time(double const seconds);
void* Computation(void *input);
#if USE_MATRIX
void Find_Distance(int tid, int num, int offset);
#endif


/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/

inline float Calculate_Distance(float *x, int i, int j, int Dim, int N_points){
	float temp=0;
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
int *label, *c_id, *local_id, *info;
float *x, *sum_dis, *min_c, *local_min;
pthread_t *pthreads;
pthread_mutex_t mutex1;
pthread_mutex_t mutex2;
pthread_barrier_t barrier1;
pthread_barrier_t barrier2;
#if USE_MATRIX
float **distance_m;
#endif


int main(int argc, char** argv){
	int i,j,num,offset, ier;
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


	for(i=0; i<N_c; i++){
		c_id[i]=i;
		min_c[i]=SUM_MAX;		
	}
	num=N_points/N_thread;
	flag=1;
	offset=0;
	//creat pthread here
	pthread_barrier_init(&barrier1, NULL, N_thread);
	pthread_barrier_init(&barrier2, NULL, N_thread);
	pthread_mutex_init(&mutex1, NULL);
	pthread_mutex_init(&mutex2, NULL);
	for(i=0; i<N_thread; i++){
		info[3*i]=i;
		info[1+3*i]=(i<(N_points-num*N_thread))?num+1:num;
		info[2+3*i]=offset;
		offset+=info[1+3*i];
		ier=pthread_create(&pthreads[i], NULL, Computation, (void*)&info[3*i]);
	//	printf("%d at %d creating\n",ier,i);
	}
	//end pthread here
	for(i=0; i<N_thread; i++){
		ier=pthread_join(pthreads[i],NULL);
	//	printf("%d at %d joining\n",ier,i);
	}
	pthread_barrier_destroy(&barrier1);
	pthread_barrier_destroy(&barrier2);
	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
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
	x=(float*)malloc(N_points*Dim*sizeof(float));
	sum_dis=(float*)malloc(N_points*sizeof(float));
	min_c=(float*)malloc(N_c*sizeof(float));
	local_min=(float*)malloc(N_thread*N_c*sizeof(float));
	label=(int*)malloc(N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
	local_id=(int*)malloc(N_thread*N_c*sizeof(int));
	info=(int*)malloc(N_thread*3*sizeof(int));
	pthreads=(pthread_t*)malloc(N_thread*N_c*sizeof(pthread_t));
	
#if USE_MATRIX
	int i;
	distance_m=(float**)malloc(N_points*sizeof(float*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(float*)malloc(N_points*sizeof(float));
	}
#endif

}
#if USE_MATRIX

void Find_Distance(int tid, int num, int offset){
	int i, j, k;
//	#pragma omp for
	for(i=offset; i<offset+num; i++){
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

void Label_Point(int tid, int num, int offset){
	int i, j, k;
	float sum;
	float temp;
//		printf("called label point:\n");
//	#pragma omp for
	for(i=offset; i<offset+num; i++){
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

int Find_Medroid(int tid, int num, int offset){
	int i,j,k;
	int temp,flag;
//	printf("called find medroid\n");
	flag=1;
	for(i=0; i<N_c; i++){
		local_min[i+tid*N_c]=min_c[i];
		local_id[i+tid*N_c]=c_id[i];
	}
//	#pragma omp for
	for(i=offset; i<offset+num; i++){
		sum_dis[i]=0;
	}
//	#pragma omp barrier
	pthread_barrier_wait(&barrier1);

//	#pragma omp for
	for(i=offset; i<offset+num; i++){
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
		if(sum_dis[i]<local_min[label[i]+tid*N_c]){
			local_id[label[i]+tid*N_c]=i;
			local_min[label[i]+tid*N_c]=sum_dis[i];
		}
	}
//	#pragma omp barrier
	pthread_barrier_wait(&barrier2);
//	#pragma omp critical
	pthread_mutex_lock( &mutex1 ); 
	for(i=0; i<N_c; i++){
		if(local_min[i+tid*N_c]<min_c[i]){
			c_id[i]=local_id[i+tid*N_c];
			min_c[i]=local_min[i+tid*N_c];
			flag=0;
		}
	}
	pthread_mutex_unlock( &mutex1 ); 	
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
	free(info);
	free(pthreads);
#if USE_MATRIX
	int i;
	for(i=0; i<N_points; i++){
		free(distance_m[i]);
	}
	free(distance_m);
#endif
}

void* Computation(void* input){
		int tid,num,offset;
		int i=0, temp=0;
		tid = *((int*)input);
		num = *((int*)input+1);
		offset = *((int*)input+2);
		printf("tid=%d, num=%d, offset=%d\n",tid,num, offset);

#if USE_MATRIX
		Find_Distance(tid, num, offset);
		pthread_barrier_wait(&barrier2);
#endif
 		Label_Point(tid, num, offset);
//		printf("%d %d at %d\n",flag, i<N_IT, tid);
		while(flag && (i<N_IT)){
			temp=Find_Medroid(tid, num, offset);

//			#pragma omp critical
//			{
			pthread_mutex_lock( &mutex2 ); 
				flag&=temp;
			pthread_mutex_unlock( &mutex2 ); 
//			}
//			#pragma omp barrier
			pthread_barrier_wait(&barrier1);
			if(tid==0){
				flag=!flag;
			}
			Label_Point(tid, num, tid);
			i++;
		}
		printf("break at it=%d at %d\n",i-1, tid);
		pthread_exit(0); 
}