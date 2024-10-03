#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <omp.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define N_IT 20
#define SUM_MAX 1e5

void Create_Memory();
char Load_File(char *str);
void Free_Memory();
int Find_Medroid(double *loacal_min, int *local_id);
void Label_Point();
void Save_Result();
void Find_Distance();

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
double *x, *sum_dis, *min_c, *distance_m;

int main(int argc, char** argv){
	int i,j,num;
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
	omp_set_dynamic(0);
	omp_set_num_threads(N_thread);

	for(i=0; i<N_c; i++){
		c_id[i]=i;
		
	}
	num=N_points/N_thread;
	#pragma omp parallel
	{
		int tid=omp_get_thread_num();
		double *local_min;
		int *local_id;
		int N_local=(tid<N_points-num*N_thread)?num+1:num;
		printf("tid=%d, N_local=%d\n",tid, N_local);

		local_min=(double*)malloc(N_c*sizeof(double));
		local_id=(int*)malloc(N_c*sizeof(double));

		Find_Distance();
		Label_Point();
		for(i=0; i<N_IT; i++){
			if(Find_Medroid(local_min,local_id)==0){
				printf("Converged!!\n");
				break;
			}
			Label_Point();
		}
		free(local_id);
		free(local_min);
	}
	
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
	for(j=0; j<Dim; j++){
		for(i=0; i<N_points; i++){
			printf("%lf ",x[i+j*N_points]);
		}
		printf("\n");
	}
		
	return 0;
}

void Find_Distance(){
	int i, j, k;
	#pragma omp for
	for(i=0;i<N_points; i++){
		for(j=0; j<N_points; j++){
			distance_m[i+j*N_points]=0;
			for(k=0; k<Dim; k++){
				distance_m[i+j*N_points]+=(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
			}
//			printf("%lf ", distance_m[i+j*N_points]);
		}
//		printf("\n");
	}

}

void Create_Memory(){
	x=(double*)malloc(N_points*Dim*sizeof(double));
	sum_dis=(double*)malloc(N_points*sizeof(double));
	min_c=(double*)malloc(N_c*sizeof(double));
	distance_m=(double*)malloc(N_points*N_points*sizeof(double));
	label=(int*)malloc(N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
	
}

void Label_Point(){
	int i, j, k;
	double sum;
	#pragma omp for
	for(i=0; i<N_points; i++){
		sum=SUM_MAX;
		label[i]=-1;
		for(j=0; j<N_c; j++){
			if(distance_m[i+c_id[j]*N_points]<sum){
				label[i]=j;
				sum=distance_m[i+c_id[j]*N_points];
			}
		}
	}
}

int Find_Medroid(double *local_min, int *local_id){
	int flag=0,temp=0;
	int i,j,k;
	#pragma omp barrier
	for(i=0; i<N_c; i++){
		local_min[i]=min_c[i];
	}
	
	#pragma omp for private(temp)
	for(i=0; i<N_points; i++){
		sum_dis[i]=0;
		for(j=0; j<N_points; j++){
			if(label[i]==label[j]){
				sum_dis[i]+=distance_m[i+j*N_points];
			}
		}
		if(sum_dis[i]<local_min[label[i]]){
			local_id[label[i]]=i;
			local_min[label[i]]=sum_dis[i];
			temp++;
		}
	}
	#pragma omp barrier
	#pragma omp critical
	{
		flag+=temp;
		for(i=0; i<N_c; i++){
			if(local_min[i]<min_c[i]){
				c_id[i]=local_id[i];
				min_c[i]=local_min[i];
			}
		}
	}
	return flag;
}

void Save_Result(){
	FILE *out;
	int i;
	out = fopen("label.txt","w");
	for(i=0; i<N_points; i++){
		fprintf(out, "%d\n",label[i]);
	}
	fclose(out);
}


void Free_Memory(){
	free(x);
	free(sum_dis);
	free(min_c);
	free(label);
	free(c_id);
}
