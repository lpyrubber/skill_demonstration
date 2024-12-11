#include "km_cuda.h"


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

static void print_time(double const seconds);
char Load_File(char *str);
void Calculate_2N();
void Initialize();



int NC, NC2, TPB, BPG, BPG2, BPGC, NP, NP2, Dim, NT, NB;


int *h_label, *h_id, *h_id_old, *h_flag;
double *h_x;

int *d_label, *d_id, *d_id_old, *d_nlist, *d_prefix, *d_flag;
double *d_x, *d_min, *d_sum;

int main(int argc, char** argv){
	double st,et;
	if(argc<4) {
		printf("Not enough of arguemt\n");
		return 1;
	}
	if(argc>4) {
		printf("Too many arguements\n");
		return 2;
	}
	
	NC=atoi(argv[2]);
	TPB= atoi(argv[3]);
	if(Load_File(argv[1])){
		return 3;
	}

	st=monotonic_seconds();

	Initialize();
    Send_To_Device();
    GPU_Compute();
    Send_To_Host();

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
	if(fscanf(in,"%d %d\n",&NP, &Dim )!=2){
		printf("can't get Number of nodes and relative dimension\n");
		return 2;
	}
    Calculate_2N();
	Allocate_Memory();
	for(i=0; i<NP; i++){
		for(j=0; j<Dim-1; j++){
			if(fscanf(in, "%f ", h_x+(j+Dim*i))<1){
				printf("can't get data\n");
				return 3;
			}
		}
		if(fscanf(in, "%f\n", h_x+((Dim-1)+Dim*i))<1){
			printf("can't get data\n");
			return 3;
		}
	}
		
	return 0;
}


void Calculate_2N(){
    NC2=0;
    NP2=1;
    while(NP2<NP){
        NP2<<=1;
    }
    while(NC2<NC){
        NC2+=TPB;
    }
    BPG=(int)((NP+TPB-1)/TPB);
    BPG2=(int)((NP2+TPB-1)/TPB);
    BPGC=(int)(NC2/TPB);
    printf("NP=%d, NP2=%d, NC=%d, NC2=%d, Dim=%d\n",NP, NP2, NC, NC2,Dim);
    printf("TPB=%d, BPG=%d, BPG2=%d, BPGC=%d\n",TPB, BPG, BPG2, BPGC);
}

void Initialize(){
    int i;
    for(i=0; i<NC; i++){
		h_id[i]=i;
		h_id_old[i]=i;
	}
}

void Save_Result(){
	FILE *out;
	int i,j;
	out = fopen("clusters.txt","w");
	for(i=0; i<NP; i++){
		fprintf(out, "%d\n",h_label[i]);
	}
	fclose(out);
	out = fopen("centroids.txt","w");
	for(i=0; i<NC; i++){
		for(j=0; j<Dim-1; j++){
			fprintf(out, "%lf ",h_x[j+h_id[i]*Dim]);
		}
		fprintf(out,"%lf\n",h_x[(Dim-1)+h_id[i]*Dim]);
	}
	fclose(out);
}
