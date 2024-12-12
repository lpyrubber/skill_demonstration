#include "km_cuda.h"

static void print_time(double const seconds);
char Load_File(char *str);
void Calculate_2N();
void Initialize();

static inline double monotonic_seconds(){
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




int NC, NC2, TPB, BPG, BPG2, BPGC, BPGR, NP, NP2, Dim, Nij;
int *h_label, *h_id, *h_id_old, *h_flag;
double *h_x;
int *d_label, *d_clist, *d_id, *d_id_old, *d_nlist, *d_prefix, *d_flag;
double *d_x, *d_min, *d_sum;

int main(int argc, char** argv){
	double st,et;
	if(argc<5) {
		printf("Not enough of arguemt\n");
		return 1;
	}
	if(argc>5) {
		printf("Too many arguements\n");
		return 2;
	}
	
	NC=atoi(argv[2]);
	BPG=atoi(argv[3]);
    TPB=atoi(argv[4]);
    
	if(Load_File(argv[1])){
		return 3;
	}

	st=monotonic_seconds();
	Initialize();
    Send_To_Device();
    GPU_Compute();
    Send_To_Host();
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
			if(fscanf(in, "%lf ", h_x+(j+Dim*i))<1){
				printf("can't get data\n");
				return 3;
			}
		}
		if(fscanf(in, "%lf\n", h_x+((Dim-1)+Dim*i))<1){
			printf("can't get data\n");
			return 3;
		}
	}
		
	return 0;
}


void Calculate_2N(){
    if(BPG<0){
        BPG=(int)((NP+TPB-1)/TPB);
        Nij=1;
    }else{
        Nij=(int)((NP+BPG*TPB-1)/(BPG*TPB));
    }
    
    //fix for sorting, prefix sum and reduction
    NC2=0;
    NP2=1;
    while(NP2<NP){
        NP2<<=1;
    }
    while(NC2<NC){
        NC2+=DTPB;
    }
    BPG2=(int)((NP2+DTPB-1)/DTPB);
    BPGC=(int)(NC2/DTPB);
    BPGR=(int)((NP+DTPB-1)/DTPB);
    //printf("NP=%d, NP2=%d, NC=%d, NC2=%d, Nij=%d, Dim=%d\n",NP, NP2, NC, NC2, Nij,Dim);
    //printf("TPB=%d, DTPB=%d, BPG=%d, BPG2=%d, BPGC=%d BPGR=%d\n",TPB, DTPB, BPG, BPG2, BPGC, BPGR);
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
	out = fopen("medoids.txt","w");
	for(i=0; i<NC; i++){
		for(j=0; j<Dim-1; j++){
			fprintf(out, "%e ",h_x[j+h_id[i]*Dim]);
		}
		fprintf(out,"%e\n",h_x[(Dim-1)+h_id[i]*Dim]);
	}
	fclose(out);
}
