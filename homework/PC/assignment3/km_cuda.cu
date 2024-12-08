#include <stdlib.h>
#include <stdio.h>
/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define N_IT 20
#define SUM_MAX 1e14
#define USE_MATRIX 1

//c function
void Allocate_Memory();
void Free_Memory();


char Load_File(char *str);
char Judge();
void Calculate_2N();
void GPU_Compute();
void Find_Medroid();
void Label_Point();
void Save_Result();


void Send_To_Device();
void Send_To_Host();

void Prefix_sum(int *d_a, int *d_b);
void Bitonic_Sort(int *d_a, char flag);

static void print_time(double const seconds);


//cuda function
__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total);
__global__ void Bitonic_Sort_Step(int *d_a, int j, int k, int N, char flag);
__global__ void Prefix_Down_Sweep(int *d_a, int *d_b);
__global__ void Offest_Between_Block(int *d_b, int N);
__global__ void Add_Offset(int *d_a, int *d_b);

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


int NC, NC2, TPB, BPG, NP, NP2, Dim, NT, NB;


int *h_label, *h_id;
float *h_x;

int *d_label, *d_id, *d_nlist, *d_prefix;
float *d_x, *d_min, *d_sum;

int main(int argc, char** argv){
	int i;
	double st,et;
	if(argc<4) {
		printf("Not enough of arguemt\n");
		return 1;
	}
	if(argc>4) {
		printf("Too many arguements\n");
		return 2;
	}
	
	NC=atoi(argv[2]);void Calculate_2N();
	NT= atoi(argv[3]);
	if(Load_File(argv[1])){
		return 3;
	}

	st=monotonic_seconds();

	for(i=0; i<NC; i++){
		h_id[i]=i;
	}
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
    NC2=1;
    NP2=0;
    while(NP2<NP){
        NP2<<=1;
    }
    while(NC2<NC){
        NC2+=TPB;
    }
}


void Allocate_Memory(){
    cudaError_t Error;
    size_t size=(NP2+1)*sizeof(int);
    
    
    
    h_x=(float*)malloc(NP*Dim*sizeof(double));
    h_label=(int*)malloc(NP*sizeof(int));
    h_id=(int*)malloc(NC*sizeof(int));


    Error=cudaMalloc((void**) &d_x, size);
    printf("CUDA Error (malloc d_a) = %s\n", cudaGetErrorString(Error));
    size=(BPG+1)*sizeof(int);
    Error=cudaMalloc((void**) &d_id, size);
    printf("CUDA Error (malloc d_b) = %s\n", cudaGetErrorString(Error));

}
void Send_To_Device(){
    cudaError_t Error;
    size_t size=(NC2+1)*sizeof(int);
    Error=cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_a -> d_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void GPU_Compute(){
    char flag=1;
    int it=0;
    while((it<N_IT)&&flag){
        Label_Point();
        Find_Medroid();
        flag=Judge();
        it++;
    }
    Bitonic_Sort(d_label,1);
}

void Send_To_Host(){
    cudaError_t Error;
    size_t size=(NC2+1)*sizeof(int);
    Error=cudaMemcpy(h_x, d_x, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_a -> h_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void Free_Memory(){
    cudaError_t Error;
    free(h_x);
    Error=cudaFree(d_x);
    printf("CUDA Error ( free d_a )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_id);
    printf("CUDA Error ( free d_b )=%s\n", cudaGetErrorString(Error));

}

char Judge(){
    return 0;
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

void Prefix_sum(int *d_a, int *d_b){
    printf("TPB=%d, BPG=%d, N2=%d, N_total=%d\n",TPB, BPG, NC2, NP);
    Prefix_Up_Sweep<<<BPG,TPB>>>(d_a, NC2, NP);
    Prefix_Down_Sweep<<<BPG,TPB>>>(d_a, d_b);
    if(BPG>1){
        Offest_Between_Block<<<1,1>>>(d_b, BPG);
        Add_Offset<<<BPG, TPB>>>(d_a, d_b);
    }
}

void Bitonic_Sort(int* d_a, char flag){
    int j,k;

    printf("TPB=%d, BPG=%d\n",TPB, BPG);
    
    for(k=2; k<=NP2; k<<=1){
        for(j=k>>1; j>0; j>>=1){
            Bitonic_Sort_Step<<<BPG,TPB>>>(d_a, j, k, NP2, flag);
        }
    }
}

void Find_Medroid(){

}

void Label_Point(){
    
}


__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
	int I = threadIdx.x;
    int l;
    if(i==N-1){
        d_a[i+1]=N_total;
    }
    __syncthreads();
    for(int stride=2; stride<=blockDim.x; stride<<=1){
        l=stride>>1;
        if((I&(stride-1))==0){
            d_a[i+stride-1]+=d_a[i+stride-1-l];
        }
        __syncthreads();
    }
    

}
__global__ void Prefix_Down_Sweep(int *d_a, int *d_b){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
	int I = threadIdx.x;
    int N2 = blockDim.x*(blockIdx.x+1);
    if(I==blockDim.x-1){
        d_b[blockIdx.x+1]+=d_a[i];
        d_a[i]=0;
    }
 
    __syncthreads();
    for(int stride=blockDim.x; stride>1; stride>>=1){
        int l=stride>>1;
        if((I&(stride-1))==0){
            int t=d_a[N2-1-I-l];
            d_a[N2-1-I-l]=d_a[N2-1-I];
            d_a[N2-1-I]+=t;
        }
        __syncthreads();
    }
}

__global__ void Offest_Between_Block(int *d_b, int N){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    //serial prefix scan
    if(i==0){
        d_b[0]=0;
        for(int i1=0; i1<N; i1++){
            d_b[i1+1]+=d_b[i1];
        }
    }
}

__global__ void Add_Offset(int *d_a, int *d_b){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    d_a[i]+=d_b[blockIdx.x];
}

__global__ void Bitonic_Sort_Step(int *d_a, int j, int k, int N, char flag){
    int i, ixj,m;
    i = threadIdx.x + blockDim.x * blockIdx.x;
    if(i<N){
        m=i&k;
        if(flag){
            i+=N;
        }
        ixj=i^j;
        
        if((ixj>i)){
            if(m==0){
                if(d_a[i]>d_a[ixj]){
                    int temp=d_a[i];
                    d_a[i]=d_a[ixj];
                    d_a[ixj]=temp;
                    int i2=(i+N)&((N<<1)-1);
                    int ixj2=(ixj+N)&((N<<1)-1);
                    temp=d_a[i2];
                    d_a[i2]=d_a[ixj2];
                    d_a[ixj2]=temp;
                }
            }else{
                if(d_a[i]<d_a[ixj]){
                    int temp=d_a[i];
                    d_a[i]=d_a[ixj];
                    d_a[ixj]=temp;
                    int i2=(i+N)&((N<<1)-1);
                    int ixj2=(ixj+N)&((N<<1)-1);
                    temp=d_a[i2];
                    d_a[i2]=d_a[ixj2];
                    d_a[ixj2]=temp;
                }
            }
        }
    }
}