#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <math.h>
#include <vector>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define N_IT 1
#define SUM_MAX 1e14
#define USE_MATRIX 0

int Partition(int *array, int low, int high);
void Quicksort(int *array, int low, int high);

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
		temp+=100*(x[i+k*N_points]-x[j+k*N_points])*(x[i+k*N_points]-x[j+k*N_points]);
	}
	return sqrt(temp);
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

void Swap(int* array, int i, int j, int N_points){
    int temp=array[i];
	int offset=(i<N_points)?N_points:-N_points;
    array[i]=array[j];
    array[j]=temp;
	temp=array[i+offset];
	array[i+offset]=array[j+offset];
    array[j+offset]=temp;	
}

int N_c, N_thread, N_points, Dim;
int *label, *c_id, *n_list, *prefix;
double *x, *sum_dis, *min_c;

#if USE_MATRIX
double **distance_m;
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
	std::vector<std::vector<int> > c_list(N_c);
//	printf("N_c =%d, N_thread = %d\n", N_c, N_thread);
	//intialization
	st=monotonic_seconds();
#if USE_MATRIX
		Find_Distance();
#endif
	for(i=0; i<N_c; i++){
		c_id[i]=i;
		min_c[i]=SUM_MAX;
	}

	i=0;
	flag=1;
	while(flag&&i<N_IT){
		double st1=monotonic_seconds();
		Label_Point();
		st1=monotonic_seconds()-st1;
		printf("%d:label_point:%lf\n",i,st1);
		double st2=monotonic_seconds();
		flag=Find_Medroid();
		st2=monotonic_seconds()-st2;
		printf("%d:find_medroid:%lf\n",i,st2);
		i++;
	}
	Quicksort(label, N_points, 2*N_points-1);
//	printf("break at %d\n",i-1);
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
	if(fscanf(in,"%d %d\n",&N_points, &Dim )!=2){
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
	prefix=(int*)malloc((N_c+1)*sizeof(int));
	n_list=(int*)malloc(N_c*sizeof(int));
#if USE_MATRIX
	int i;
	distance_m=(double**)malloc(N_points*sizeof(double*));
	for(i=0;i<N_points;i++){
		distance_m[i]=(double*)malloc((N_points-i)*sizeof(double));
	}
#endif	
}
#if USE_MATRIX

void Find_Distance(){
	int i, j, k;
	#pragma omp for
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
void Label_Point(){
//	printf("start label\n");
	int i, j, k, ip, jp;
	double sum,temp;
	
	for(i=0; i<N_c;i++){
		n_list[i]=0;
	}
	for(i=0; i<N_points; i++){
		sum=SUM_MAX;
		label[i]=-1;
		label[i+N_points]=i;
		for(j=0; j<N_c; j++){
#if USE_MATRIX
			if(c_id[j]>i){
				ip=i;
				jp=c_id[j];
			}else{
				ip=c_id[j];
				jp=i;
			}
			if(distance_m[ip][jp-ip]<sum){
				label[i]=j;
				sum=distance_m[ip][jp-ip];
			}

#else			
			temp=Calculate_Distance(x,i,c_id[j],Dim,N_points);
			if(i==796){
				printf("%d: C=%d with dis=%lf, label=%d, and sum=%lf\n",i,c_id[j],temp, label[i], sum);
			}
			if(temp<sum){
				label[i]=j;
				sum=temp;
			}
#endif
		}
		n_list[label[i]]++;	
	}
	Quicksort(label, 0, N_points-1);
//	printf("Finish quicksort\n");
//	for(i=0;i<N_points;i++){
//		printf("%d, %d\n",label[i],label[i+N_points]);
//	}
	prefix[0]=0;
	for(i=0; i<N_c;i++){
		prefix[i+1]=prefix[i]+n_list[i];
	}
/*	
	for(i=0;i<N_c+1;i++){
		printf("prefix[%d]=%d\n", i, prefix[i]);
	}
	k=1;
	for(i=0;i<N_points;i++){
		if(prefix[k]==i){
			printf("-------------------------\n");
			k++;
		}
		printf("%d: %d, %d\n",i, label[i],label[i+N_points]);
	}
*/
}

int Find_Medroid(){
	int flag=0;
	int i,j,k, ip, jp, id_new, im, jm;
	double temp=0, sum, min;
	FILE *in;
	in = fopen("serial_sum.csv","w");
	for(i=0; i<N_c; i++){
		min=SUM_MAX;
		for(j=prefix[i]; j<prefix[i+1]; j++){
			ip=label[j+N_points];
			sum=0;
			for(k=prefix[i];k<prefix[i+1];k++){
				jp=label[k+N_points];
			
#if USE_MATRIX
				im=(jp>ip)?ip:jp;
				jm=(jp>ip)?jp:ip;
				sum+=distance_m[im][jm-im]/n_list[i];
#else
				temp=Calculate_Distance(x,ip,jp,Dim,N_points)/n_list[i];		
				sum+=temp;
#endif						

			}
			fprintf(in, "%d,%f,%d\n",ip,sum,label[ip]);
			if(sum<min){
				id_new=ip;
				min=sum;
			}
		}
		if(id_new!=c_id[i]){
			flag++;
			c_id[i]=id_new;
		}
	}
	fclose(in);
	return (flag>0)?1:0;
}

void Save_Result(){
	FILE *out;
	int i,j;
	out = fopen("clusters.txt","w");
	for(i=0; i<N_points; i++){
		fprintf(out, "%d\n",label[i]);
	}
	fclose(out);
	out = fopen("medoids.txt","w");
	for(i=0; i<N_c; i++){
		fprintf(out, "%d\n", c_id[i]);
/*		for(j=0; j<Dim-1; j++){
			fprintf(out, "%lf ",x[c_id[i]+j*N_points]);
		}
		fprintf(out,"%lf\n",x[c_id[i]+(Dim-1)*N_points]);
*/	}
	fclose(out);
}

void Quicksort(int *array, int low, int high){
    int pi,tag;
    if (low < high){
        pi = Partition(array, low, high);
        Quicksort(array, low, pi-1);
        Quicksort(array, pi+1, high);
    }
}

int Partition(int *array, int low, int high){

	int pivot=array[high];
    int j = low;
    int i = (low-1);
    for(j=low; j<= high-1; j++){
        if(array[j]<pivot){   
            i++;
            Swap(array, i, j, N_points);
        }
    }
    Swap(array, i+1, high, N_points);
    return i+1;
}

void Free_Memory(){
	free(x);
	free(sum_dis);
	free(min_c);
	free(label);
	free(c_id);
	free(n_list);
	free(prefix);
#if USE_MATRIX
	int i;
	for(i=0; i<N_points; i++){
		free(distance_m[i]);
	}
	free(distance_m);
#endif
}
