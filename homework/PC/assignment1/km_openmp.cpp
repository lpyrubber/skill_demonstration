#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <omp.h>
#include <math.h>
#include <vector>
/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define USE_MATRIX 1
#define N_IT 20
#define SUM_MAX 1e14

void Create_Memory();
char Load_File(char *str);
void Free_Memory();
void Find_Medroid(int tid);
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

int N_c, N_thread, N_points, Dim, flag, c_id_new;
int *label, *c_id, *local_id;
float *x, *sum_dis, *min_c, *local_min;
#if USE_MATRIX
float **distance_m;
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
	std::vector<std::vector<int> > c_list(N_c);
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
		int i,j,k,ip,jp,im,jm;
		int itt=0;
		double temp=0, sum, min;
		std::vector<std::vector<int>> local_list(N_c);
		
#if USE_MATRIX
		Find_Distance();
#endif

		while(flag && (itt<N_IT)){
			#pragma omp barrier
			#pragma omp for
			for(i=0; i<N_points; i++){
				min=SUM_MAX;
				label[i]=-1;
				for(j=0; j<N_c; j++){
#if USE_MATRIX
					im=(c_id[j]>i)?i:c_id[j];
					jm=(c_id[j]>i)?c_id[j]:i;
					if(distance_m[im][jm-im]<min){
						label[i]=j;
						min=distance_m[im][jm-im];
					}

#else			
					temp=Calculate_Distance(x,i,c_id[j],Dim,N_points);
					if(temp<min){
						label[i]=j;
						min=temp;
					}
#endif
				}
			}
			
			#pragma omp for schedule(static)
			for (i = 0; i < N_points; i++) {
                local_list[label[i]].push_back(i);
            }
			#pragma omp master
			{
				for(i=0;i<N_c;i++){
					c_list[i].clear();
				}
			}
			#pragma omp barrier
			#pragma omp critical
			{
				for(i=0; i<N_c; i++){
					c_list[i].insert(c_list[i].end(), local_list[i].begin(), local_list[i].end());
					local_list[i].clear();
				}
			}
			flag=0;
			#pragma omp barrier
/*			#pragma  omp master
			{
			for(i=0; i<N_c; i++){
				printf("at cluaster %d, N_c=%d\n",i,c_list[i].size());
				for(j=0; j<c_list[i].size();j++){
					printf("%d ",c_list[i][j]);
				}
				printf("\n\n");
			}			
*/
			
			
			for(i=0; i<N_c; i++){
				min_c[i]=SUM_MAX;
				#pragma omp for
				for(j=0; j<c_list[i].size(); j++){
					ip=c_list[i][j];
					sum=0;
					for(k=0; k<c_list[i].size(); k++){
						jp=c_list[i][k];
#if USE_MATRIX
						im=(jp>ip)?ip:jp;
						jm=(jp>ip)?jp:ip;
						sum+=distance_m[im][jm-im]/c_list[i].size();
#else
						temp=Calculate_Distance(x,ip,jp,Dim,N_points)/c_list[i].size;		
						sum+=temp;
#endif				
					}
					#pragma omp critical
					{
						if(sum<min_c[i]){
							min_c[i]=sum;
							c_id_new=ip;
						}
					}
				}
				#pragma omp barrier
				#pragma omp master
				{
					if(c_id_new!=c_id[i]){
						c_id[i]=c_id_new;
						flag++;
					}
				}
				#pragma omp barrier
			}
			itt++;
		}
//		printf("%d,break at it=%d\n",itt,tid);	
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
	x=(float*)malloc(N_points*Dim*sizeof(float));
	sum_dis=(float*)malloc(N_points*sizeof(float));
	min_c=(float*)malloc(N_c*sizeof(float));
	local_min=(float*)malloc(N_thread*N_c*sizeof(float));
	label=(int*)malloc(N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
	local_id=(int*)malloc(N_thread*N_c*sizeof(int));
#if USE_MATRIX
	int i;
	distance_m=(float**)malloc(N_points*sizeof(float*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(float*)malloc((N_points-i)*sizeof(float));
	}
#endif

}
#if USE_MATRIX
void Find_Distance(){
	int i, j, k;
	#pragma omp for schedule(dynamic)
	for(i=0; i<N_points; i++){
		for(j=i; j<N_points; j++){
			distance_m[i][j-i]=0;
			for(k=0; k<Dim; k++){
				distance_m[i][j-i]+=(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
			}
			distance_m[i][j-i]=sqrtf(distance_m[i][j-i]);
		}
	}
}
#endif

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
