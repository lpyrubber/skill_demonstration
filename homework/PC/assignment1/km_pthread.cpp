#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <pthread.h>
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
int Find_Medroid( int tid, int time);
void Label_Point( int tid,   int num,   int offset);
void Save_Result();
static void print_time(double const seconds);
void* Computation(void *input);
#if USE_MATRIX
void Find_Distance( int tid,   int num,   int offset);
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

int N_c, N_thread, N_points, Dim, flag, c_new;
int *label, *c_id, *info, *c_old;
float *x, *sum_dis, *min_c;
pthread_t *pthreads;
pthread_mutex_t mutex1;
pthread_mutex_t mutex2;
pthread_barrier_t barrier1;
pthread_barrier_t barrier2;
#if USE_MATRIX
float **distance_m;
#endif
std::vector<std::vector<int> > c_list;


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
//	printf("N_c =%d, N_thread = %d\n", N_c, N_thread);
    c_list.resize(N_c);
	st=monotonic_seconds();
    

	for(i=0; i<N_c; i++){
		c_id[i]=i;
		min_c[i]=SUM_MAX;		
	}
	num=N_points/N_thread;
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
//		printf("%d at %d creating\n",ier,i);
	}
	//end pthread here
	for(i=0; i<N_thread; i++){
		ier=pthread_join(pthreads[i],NULL);
//		printf("%d at %d joining\n",ier,i);
	}
	pthread_barrier_destroy(&barrier1);
	pthread_barrier_destroy(&barrier2);
	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
	et=monotonic_seconds();
	print_time(et-st);
	Save_Result();
//    printf("finish save\n");
	Free_Memory();
//    printf("finish free\n");
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
	if(fscanf(in,"%d %d\n",&N_points, &Dim )!=2){
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
	label=(int*)malloc(N_points*sizeof(int));
	c_id=(int*)malloc(N_c*sizeof(int));
    c_old=(int*)malloc(N_c*sizeof(int));
	info=(int*)malloc(N_thread*3*sizeof(int));
	pthreads=(pthread_t*)malloc(N_thread*sizeof(pthread_t));
	
#if USE_MATRIX
	int i;
	distance_m=(float**)malloc(N_points*sizeof(float*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(float*)malloc((N_points-1)*sizeof(float));
	}
#endif

}
#if USE_MATRIX

void Find_Distance( int tid, int num, int offset){
	int i, j, k, ib, ix, iy;
	for(i=offset; i<offset+num; i++){
		ib=(i%2)?((i%2)*(N_points-1)-(i/2)):((i%2)*(N_points-1)+(i/2));
		for(j=ib; j<N_points; j++){
			distance_m[ib][j-ib]=0;
			for(k=0; k<Dim; k++){
				distance_m[ib][j-ib]+=(x[ib+k*N_points]-x[j+k*N_points])*(x[ib+k*N_points]-x[j+k*N_points]);
			}
			distance_m[ib][j-ib]=sqrtf(distance_m[ib][j-ib]);
		}
	}

}
#endif

void Label_Point( int tid,  int num,  int offset){
	int i, j, k;
    int im,jm;
	float min;
	float temp;
	std::vector<std::vector<int>> local_list(N_c);
    if(tid==0){
        for(i=0;i<N_c;i++){
            c_list[i].clear();
        }
    }
	for(i=offset; i<(offset+num); i++){
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
        local_list[label[i]].push_back(i);

	}
    pthread_barrier_wait(&barrier2);
    pthread_mutex_lock( &mutex1 );
	for(i=0; i<N_c; i++){
		c_list[i].insert(c_list[i].end(), local_list[i].begin(), local_list[i].end());
		local_list[i].clear();
	}
	pthread_mutex_unlock( &mutex1 ); 
    if(tid==0){
        for(i=0;i<N_c;i++){
            c_old[i]=c_id[i];
        }
    }

}

int Find_Medroid( int tid, int time){
	int i,j,k, ip, jp, id_new, im, jm, offset, num, ave, flag;
	double temp=0, sum, min, local;
	flag=0;
	for(i=0; i<N_c; i++){
        min_c[i]=SUM_MAX;
		min=SUM_MAX;
        ave=c_list[i].size()/N_thread;
        num=(tid<(c_list[i].size()-ave*N_thread))?ave+1:ave;
        offset=(tid<(c_list[i].size()-ave*N_thread))?(ave+1)*tid:c_list[i].size()-(N_thread-tid)*ave;
    	for(j=offset; j<offset+num; j++){
	    	ip=c_list[i][j];
		   	sum=0;
		    for(k=0; k<c_list[i].size(); k++){
			   	jp=c_list[i][k];
#if USE_MATRIX
	    		im=(jp>ip)?ip:jp;
		   		jm=(jp>ip)?jp:ip;
		    	sum+=distance_m[im][jm-im]/c_list[i].size();
#else
			   	temp=Calculate_Distance(x,ip,jp,Dim,N_points)/c_list[i].size();		
			   	sum+=temp;
#endif				
		   	}

		    if(sum<min){
			   	min=sum;
			    local=ip;
			}

        }
		pthread_mutex_lock( &mutex2 );
		if(min<min_c[i]){
			min_c[i]=min;
			c_id[i]=local;
		}
		pthread_mutex_unlock( &mutex2 ); 
	}
    pthread_barrier_wait(&barrier1);
    for(i=0; i<N_c; i++){
        if(c_old[i]!=c_id[i]){
            flag++;
        }
    }
	return flag;
}

void Save_Result(){
	FILE *out;
	int i, j;
	out = fopen("clusters.txt","w");
	for(i=0; i<N_points; i++){
		fprintf(out, "%d\n",label[i]);
	}
	fclose(out);
	out = fopen("centroids.txt","w");
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
	free(info);
    free(c_old);
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
		int itt=0, temp=0, i;
		int tid = *((int*)input);
		int num = *((int*)input+1);
		int offset = *((int*)input+2);
//		printf("tid=%d, num=%d, offset=%d\n",tid,num, offset);
        int flag=1;
#if USE_MATRIX
		Find_Distance(tid, num, offset);
		pthread_barrier_wait(&barrier2);
#endif
		while(flag && (itt<N_IT)){   
            Label_Point(tid, num, offset);
            pthread_barrier_wait(&barrier1);
			flag=Find_Medroid(tid,itt);
			itt++;
		}
//		printf("break at it=%d at %d\n",itt, tid);
		pthread_exit(0); 
}