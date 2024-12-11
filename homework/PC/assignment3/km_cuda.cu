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

#define N_IT 1
#define SUM_MAX 1e14
#define DTPB 512

//cpu function
void Allocate_Memory();
void Free_Memory();
void Send_To_Device();
void GPU_Compute();
void Send_To_Host();
void Save_Result();
void Prefix_sum(int *d_a, int *d_b);
void Bitonic_Sort(int *d_a, char flag);
void Find_Medroid();
void Label_Point();
char Judge();
static void print_time(double const seconds);
char Load_File(char *str);
void Calculate_2N();
void Initialize();



//cuda function (didn't fintout a way to split)

__global__ void Bitonic_Sort_Step(int *d_a, int j, int k, int N, char flag);
__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total);
__global__ void Prefix_Down_Sweep(int *d_a, int *d_b);
__global__ void Offest_Between_Block(int *d_b, int N);
__global__ void Add_Offset(int *d_a, int *d_b);
__global__ void Nearest_Medroid(double *x, int *label, int *clist, int *id, int *id_old, int *nlist, int Nij, int NC, int N, int N2, int Dim);
__global__ void Nlist_Initial(int *nlist, int *label, int N);
__global__ void Nlist_Reduction(int *nlist, int offset, int gap, int NC);
__global__ void Nlist_To_Prefix(int *nlist, int *prefix, int NC, int NC2, int BBPG, int BPG, int N);
__global__ void Distance_Between_Cluster(double *x, double *sum, int *label, int *prefix, int Dim, int Nij, int N, int N2);
__global__ void Min_Initial(double *min, double *sum, int *label, int *prefix, int *id, int NC, int N, int N2);
__global__ void Min_Reduction(double *min, int *id, int offset, int gap, int NC);
__global__ void Min_Further(double *min, int *id, int NC, int BBPG, int BPG);
__global__ void Update_ID(int *id, int*id_old, int NC, int *d_flag);
__device__ double Calculate_Distance(double *x, int i, int j, int Dim);


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




int NC, NC2, TPB, BPG, BPG2, BPGC, BPGR, NP, NP2, Dim, Nij, Ncj;
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
	TPB= atoi(argv[3]);
    BPG=atoi(argv[4]);
    
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
    
    Nij=(int)((NP+BPG*TPB-1)/(BPG*TPB));
    Ncj=(int)((NC+BPG*TPB-1)/(BPG*TPB));
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
    printf("NP=%d, NP2=%d, NC=%d, NC2=%d, Nij=%d, Ncj=%d, Dim=%d\n",NP, NP2, NC, NC2, Nij, Ncj,Dim);
    printf("TPB=%d, DTPB=%d, BPG=%d, BPG2=%d, BPGC=%d BPGR=%d\n",TPB, DTPB, BPG, BPG2, BPGC, BPGR);
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


void Allocate_Memory(){   
    h_flag=(int*)malloc(BPG*sizeof(int));
    h_label=(int*)malloc(NP*sizeof(int));
    h_id=(int*)malloc(NC*sizeof(int));
    h_id_old=(int*)malloc(NC*sizeof(int));
    h_x=(double*)malloc(NP*Dim*sizeof(double));
    
     
    size_t size;

    size=NP*sizeof(int);
    cudaMalloc((void**) &d_label, size);
    size=2*NP2*sizeof(int);
    cudaMalloc((void**) &d_clist, size);
    //size=BPGR*NC*sizeof(int);
    size=NP*sizeof(int);
    cudaMalloc((void**) &d_id, size);
    size=NC*BPG*sizeof(int);
    cudaMalloc((void**) &d_id_old, size);
    size=BPGR*NC*sizeof(int);
    cudaMalloc((void**) &d_nlist, size);
    size=(NC2+1)*sizeof(int);
    cudaMalloc((void**) &d_prefix, size);
    size=BPG*sizeof(int);
    cudaMalloc((void**) &d_flag, size);
    size=NP*Dim*sizeof(double);
    cudaMalloc((void**) &d_x, size);
    size=NP*sizeof(double);
    cudaMalloc((void**) &d_sum, size);
    size=BPGR*NC*sizeof(double);
    cudaMalloc((void**) &d_min, size);
  
}
void Send_To_Device(){
    size_t size;
    size=NP*Dim*sizeof(double);
    cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);
    size=NC*sizeof(int);
    cudaMemcpy(d_id, h_id, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_id_old, h_id_old, size, cudaMemcpyHostToDevice);
}
void GPU_Compute(){
    char flag=1;
    int it=0;
    while((it<N_IT)&&flag){
        Label_Point();
        Find_Medroid();
    //    flag=Judge();
        it++;
    }
}

void Send_To_Host(){
    size_t size;

    size=NC*sizeof(int);
    cudaMemcpy(h_id, d_id, size, cudaMemcpyDeviceToHost);
    size=NP*sizeof(int);
    cudaMemcpy(h_label, d_label, size, cudaMemcpyDeviceToHost);

}
void Free_Memory(){
    free(h_label);
    free(h_id);
    free(h_id_old);
    free(h_x);

    cudaFree(d_label);
    cudaFree(d_id);
    cudaFree(d_id_old);
    cudaFree(d_nlist);
    cudaFree(d_prefix);
    cudaFree(d_flag);
    cudaFree(d_x);
    cudaFree(d_min);
    cudaFree(d_sum);    

}

char Judge(){
    
    char flag=0;
    int i;
    size_t size=sizeof(int);
    //bind with reduction (TPB fix)
    Update_ID<<<BPGR,DTPB>>>(d_id, d_id_old, NC, d_flag);
    cudaMemcpy(h_flag, d_flag, size, cudaMemcpyDeviceToHost);
    for(i=0; i<BPG; i++){
        flag|=h_flag[i];
    }
    return flag;
}


void Prefix_sum(int *d_a, int *d_b){
    printf("TPB=%d, BPGC=%d, NC2=%d, N_total=%d\n",TPB, BPGC, NC2, NP);
    Prefix_Up_Sweep<<<BPGC,DTPB>>>(d_a, NC2, NP);
    Prefix_Down_Sweep<<<BPGC,DTPB>>>(d_a, d_b);
    if(BPG>1){
        Offest_Between_Block<<<1,1>>>(d_b, BPG);
        Add_Offset<<<BPGC, DTPB>>>(d_a, d_b);
    }
}

void Bitonic_Sort(int* d_a, char flag){
    int j,k;   
    for(k=2; k<=NP2; k<<=1){
        for(j=k>>1; j>0; j>>=1){
            Bitonic_Sort_Step<<<BPG2,DTPB>>>(d_a, j, k, NP2, flag);
        }
    }
}

void Label_Point(){
    int BBPG=(int)((BPGR+DTPB-1)/DTPB);
    int i, offset;

    Nearest_Medroid<<<BPG,TPB>>>(d_x, d_label, d_clist ,d_id, d_id_old, d_nlist, Nij, NC,NP, NP2,Dim);           
    //reduction for number of each cluster (fix TPB)
    Nlist_Initial<<<BPGR,DTPB>>>(d_nlist, d_label, NP);
    for(i=0; i<NC; i++){
        offset=i*BPGR;
        Nlist_Reduction<<<BBPG,DTPB>>>(d_nlist, offset, BPGR, NC);
    }
    //gather to prefix (fix TPB)
    Nlist_To_Prefix<<<BPGC,DTPB>>>(d_nlist, d_prefix, NC, NC2, BBPG, BPGR, NP);
    //prefix_sum (fix TPB)
    Prefix_sum(d_prefix, d_nlist);
    //sort (fix TPB)
    Bitonic_Sort(d_clist, 0);
}

void Find_Medroid(){
    int BBPG=(int)((BPGR+DTPB-1)/DTPB);
    int i, offset; 
    int size;
    double *h_temp;
    int *h_clist, *h_tempi;

    Distance_Between_Cluster<<<BPG, TPB>>>(d_x, d_sum, d_clist, d_prefix, Dim, Nij, NP, NP2);
    
    size=NP*sizeof(double);
    h_temp=(double*)malloc(size);
    size=2*NP2*sizeof(int);
    h_clist=(int*)malloc(size);
    h_tempi=(int*)malloc(size);

    //reduction for min (fix TPB)
    Min_Initial<<<BPGR,DTPB>>>(d_min, d_sum, d_clist, d_prefix, d_id, NC, NP, NP2);
    size=BPGR*NC*sizeof(int);
    cudaMemcpy(h_tempi, d_id, size, cudaMemcpyDeviceToHost);
    size=BPGR*NC*sizeof(double);
    cudaMemcpy(h_temp, d_min, size, cudaMemcpyDeviceToHost);
    printf("first\n");
    for(i=0; i<BPGR*NC; i++){
        printf("%d: %d %f\n",i, h_tempi[i], h_temp[i]);
        if(i%(BPGR)==BPGR-1){
            printf("=====================================\n");
        }
    }  


    for(i=0; i<NC; i++){
        offset=i*BPGR;
        Min_Reduction<<<BBPG,DTPB>>>(d_min, d_id, offset, BPGR, NC);
    }
    printf("\nsecond\n");
    size=BPGR*NC*sizeof(int);
    cudaMemcpy(h_tempi, d_id, size, cudaMemcpyDeviceToHost);
    size=BPGR*NC*sizeof(double);
    cudaMemcpy(h_temp, d_min, size, cudaMemcpyDeviceToHost);
    for(i=0; i<BPGR*NC; i++){
        printf("%d: %d %f\n",i, h_tempi[i], h_temp[i]);
        if(i%(BPGR)==BPGR-1){
            printf("=====================================\n");
        }
    }  


    Min_Further<<<BPGC,DTPB>>>(d_min, d_id, NC, BBPG, BPGR);

}




__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
	int I = threadIdx.x;
    int l;

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
    if(i<gridDim.x){
        d_b[i]=0;
    }
    __syncthreads();
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

__global__ void Nearest_Medroid(double *x,int *label, int *clist,int *id, int *id_old, int *nlist, int Nij, int NC, int N, int N2, int Dim){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = Nij*i;
    double min2,temp;
    
    for(int j1=0; j1<Nij; j1++){
        if((j1+j)<NC*gridDim.x){
            nlist[j1+j]=0;
            id_old[j1+j]=id[(j1+j)%NC];
        }
    }
    __syncthreads();
    
    for(int j1=0; j1<Nij; j1++){
        if((j1+j)<N){
            min2=SUM_MAX;
            //initial label for array
            for(int i1=0; i1<NC; i1++){
                temp=Calculate_Distance(x,j+j1,id[i1],Dim);
                if(temp<min2){
                    label[j1+j]=i1;
                    min2=temp;
                }
            }
//            atomicAdd(&nlist[blockIdx.x+label[j1+j]*gridDim.x],1);
            clist[j1+j]=label[j1+j];
            clist[j1+j+N2]=j1+j;
        }
        if(j1+j<(N2-N)){
            //initial label for redundent
            clist[j1+j+N]=NC+2;
            clist[j1+j+N+N2]=N+1;
        }
    }
}

__global__ void Nlist_Initial(int *nlist, int *label, int N){
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if(i<N){
        atomicAdd(&nlist[blockIdx.x+label[i]*gridDim.x],1);
    }
}

__global__ void Nlist_Reduction(int *nlist, int offset, int gap, int NC){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int I = threadIdx.x;
    if(i<gap){
        for(int stride = blockDim.x>>1; stride>0; stride>>=1){
            if(I<stride){
                if(i+stride<gap){
                    nlist[i+offset]+=nlist[i+stride+offset];
                }
            }
            __syncthreads();
        }
    }
}

__global__ void Nlist_To_Prefix(int *nlist, int *prefix, int NC, int NC2, int BBPG,int BPG, int N){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if(i<NC2){
        prefix[i]=0;
    }
    if(i<NC){
        for(int i1=0; i1 < BBPG; i1++){
            prefix[i]+=nlist[i*BPG+i1*BBPG];
        }
    }
    if(i==0){
        prefix[NC2]=N;
    }
}

__global__ void Distance_Between_Cluster(double *x, double *sum, int *label, int *prefix, int Dim, int Nij, int N, int N2){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = i*Nij;
    for(int j1=0; j1<Nij; j1++){
        if((j1+j)<N){
            sum[j1+j]=0;
            for(int i1=prefix[label[j1+j]]; i1<prefix[label[j1+j]+1]; i1++ ){
                sum[j1+j]+=Calculate_Distance(x,label[j1+j+N2],label[i1+N2],Dim)/(prefix[label[j1+j]+1]-prefix[label[j1+j]]);
            }
        }
    }
}

__global__ void Min_Initial(double *min, double *sum, int *label, int *prefix, int *id, int NC, int N, int N2){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int I = threadIdx.x;
    int up,down;
    __shared__ double temp[1024];
    __shared__ int tag[1024];

    if(i<NC*gridDim.x){
        id[i]=N+2;
    }
    __syncthreads();
    up=0;
    down=-1;
    for(int i1=0; i1<NC; i1++){
        if((i-I)>=prefix[i1]){
            down++;
        }
        if((i-I+blockDim.x)>prefix[i1]){
            up++;
        }
    }

    //no atomic min function for double, therefore doing reduction inside a block for each cluster; 

    for(int i1=down; i1<up; i1++){
        if(i<N){
            if(label[i]==i1){
                tag[I]=label[i+N2];
                temp[I]=sum[i];
            }else{
                tag[I]=N+2;
                temp[I]=1e16;
            }         
            __syncthreads();
            for(int stride = blockDim.x>>1; stride>0; stride>>=1){
                if(I<stride){
                    if(i+stride<N){
                        if((tag[I]!=N+2)&&(tag[I+stride]!=N+2)){
                            if(temp[I+stride]<temp[I]){
                                temp[I]=temp[I+stride];
                                tag[I]=tag[I+stride];
                            }
                        }else if(tag[I+stride]!=N+2){
                            temp[I]=temp[I+stride];
                            tag[I]=tag[I+stride];
                        } 
                    }
                }
                __syncthreads();
            }

            if(I==0){
                min[blockIdx.x+i1*gridDim.x]=temp[0];
                id[blockIdx.x+i1*gridDim.x]=tag[0];
            }
          
        }
    }

}

__global__ void Min_Reduction(double *min, int *id, int offset, int gap, int NC){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int I = threadIdx.x;
    if(i<gap){
        for(int stride = blockDim.x>>1; stride>0; stride>>=1){
            if(I<stride){
                if(i+stride<gap){
                    if(min[i+offset+stride]<min[i+offset]){
                        min[i+offset]=min[i+offset+stride];
                        id[i+offset]=id[i+offset+stride];
                    }
                }
            }
            __syncthreads();
        }
    }
}

__global__ void Min_Further(double *min, int *id, int NC, int BBPG,int BPG){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if(i<NC){
        for(int i1=1; i1 < BBPG; i1++){
            if(min[i*BPG+i1*BBPG]<min[i]){
                min[i]=min[i*BPG+i1*BBPG];
                id[i]=id[i*BPG+i1*BBPG];
            }
        }
        __syncthreads();
        min[i]=min[i*BPG];
        id[i]=id[i*BPG];
    }
 }

 __global__ void Update_ID(int *id, int*id_old, int NC, int *d_flag){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if(i<gridDim.x){
        d_flag[i]=0;
    }
    __syncthreads();
    if(i<NC){
        if(id_old[i]!=id[i]){
            atomicAdd(&d_flag[blockIdx.x],1);
            id_old[i]=id[i];
        }
    }
 }

__device__ double Calculate_Distance(double *x, int i, int j, int Dim){
    double temp=0;
    int k=0;
    for(k=0; k<Dim; k++){
        temp+=(x[k+i*Dim]-x[k+j*Dim])*(x[k+i*Dim]-x[k+j*Dim]);
	}
	return sqrtf(temp);
}