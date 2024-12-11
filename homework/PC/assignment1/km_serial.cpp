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

#define N_IT 20
#define SUM_MAX 1e14
#define USE_MATRIX 0

void Create_Memory();
char Load_File(char *str);
void Free_Memory();
int Find_Medroid(std::vector<std::vector<int> >& c_list);
void Label_Point(std::vector<std::vector<int> >& c_list );
void Save_Result();
static void print_time(double const seconds);
#if USE_MATRIX
void Find_Distance();
#endif

inline double Calculate_Distance(double *x, int i, int j, int Dim, int N_points){
	double temp=0;
	int k;
	for(k=0; k<Dim; k++){
		temp+=(x[k+i*Dim]-x[k+j*Dim])*(x[k+i*Dim]-x[k+j*Dim]);
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
int *label, *c_id, *n_list, *c_id_old;
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
	std::vector<std::vector<int> > c_list(N_c);
//	printf("N_c =%d, N_thread = %d\n", N_c, N_thread);
	//intialization
	st=monotonic_seconds();
#if USE_MATRIX
		Find_Distance();
#endif
	for(i=0; i<N_c; i++){
		c_id[i]=i;
		c_id_old[i]=i;
		min_c[i]=SUM_MAX;
	}
	Label_Point(c_list);
	i=0;
	flag=1;
	while(flag&&i<1){
		double st1=monotonic_seconds();
		flag=Find_Medroid(c_list);
		st1=monotonic_seconds()-st1;
		printf("%d:find_medroid:%lf\n",i,st1);
		double st2=monotonic_seconds();
		Label_Point(c_list);
		st2=monotonic_seconds()-st2;
		printf("%d:label_point:%lf\n",i,st2);
		i++;
	}
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
			if(fscanf(in, "%lf ", x+(j+i*Dim))<1){
				printf("can't get data\n");
				return 3;
			}
		}
		if(fscanf(in, "%lf\n", x+(Dim-1+i*Dim))<1){
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
	c_id_old=(int*)malloc(N_c*sizeof(int));
	n_list=(int*)malloc(N_c*sizeof(int));
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
	#pragma omp for
	for(i=0; i<N_points; i++){
		for(j=i; j<N_points; j++){
			distance_m[i][j-i]=0;
			for(k=0; k<Dim; k++){
				distance_m[i][j-i]+=(x[k+i*Dim]-x[k+j*Dim])*(x[k+i*Dim]-x[k+j*Dim]);
			}
			distance_m[i][j-i]=sqrtf(distance_m[i][j-i]);
		}
	}
}
#endif
void Label_Point(std::vector<std::vector<int> >& c_list){
	int i, j, k, ip, jp;
	double sum,temp;
	for(i=0; i<N_c;i++){
		c_list[i].clear();
		n_list[i]=0;
	}
	for(i=0; i<N_points; i++){
		sum=SUM_MAX;
		label[i]=-1;
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
			if(temp<sum){
				label[i]=j;
				sum=temp;
			}
#endif
		}
		n_list[label[i]]++;	
		c_list[label[i]].push_back(i);
	}

}

int Find_Medroid(std::vector<std::vector<int> >& c_list){
	int flag=0;
	int i,j,k, ip, jp, id_new, im, jm;
	double temp=0, sum, min;
	FILE *in;
	in = fopen("serial_sum.csv","w");
	for(i=0; i<N_c; i++){
		min_c[i]=SUM_MAX;
	}
	for(i=0; i<N_points; i++){
		ip=i;
		sum=0;
		for(k=0;k<c_list[label[i]].size();k++){
			jp=c_list[label[i]][k];
			
#if USE_MATRIX
			im=(jp>ip)?ip:jp;
			jm=(jp>ip)?jp:ip;
			sum+=distance_m[im][jm-im]/c_list[label[i]].size();
#else
			temp=Calculate_Distance(x,ip,jp,Dim,N_points)/c_list[label[i]].size();		
			sum+=temp;
#endif						
		}
		fprintf(in,"%d,%f,%d\n",i,sum,label[i]);
		if(sum<min_c[label[i]]){
			c_id[label[i]]=ip;
			min_c[label[i]]=sum;
		}
	}
	fclose(in);
	for(i=0; i<N_c; i++){
		if(c_id_old[i]!=c_id[i]){
			flag++;
			c_id_old[i]=c_id[i];
		}
	}
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
	out = fopen("centroids.txt","w");
	for(i=0; i<N_c; i++){
		for(j=0; j<Dim-1; j++){
			fprintf(out, "%lf ",x[j+c_id[i]*Dim]);
		}
		fprintf(out,"%lf\n",x[(Dim-1)+c_id[i]*Dim]);
	}
	fclose(out);
}


void Free_Memory(){
	free(x);
	free(sum_dis);
	free(min_c);
	free(label);
	free(c_id);
	free(c_id_old);
	free(n_list);
#if USE_MATRIX
	int i;
	for(i=0; i<N_points; i++){
		free(distance_m[i]);
	}
	free(distance_m);
#endif
}
