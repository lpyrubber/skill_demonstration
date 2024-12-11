#include "km_cuda.h"


void Prefix_sum(int *d_a, int *d_b);
void Bitonic_Sort(int *d_a, char flag);
void Find_Medroid();
void Label_Point();
char Judge();



//cuda function
__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total);
__global__ void Bitonic_Sort_Step(int *d_a, int j, int k, int N, char flag);
__global__ void Prefix_Down_Sweep(int *d_a, int *d_b);
__global__ void Offest_Between_Block(int *d_b, int N);
__global__ void Add_Offset(int *d_a, int *d_b);
__global__ void Nearest_Medroid(float *x, float *min, int *label,int *id, int *nlist, int NC, int N, int N2, int Dim);
__global__ void Nlist_Reduction(int *nlist, int offset, int gap, int NC);
__global__ void Nlist_To_Prefix(int *nlist, int *prefix, int NC, int NC2, int BBPG, int BPG, int N);
__global__ void Distance_Between_Cluster(float *x, float *min, float *sum, int *label, int *prefix, int *id, int Dim, int NC, int N, int N2);
__global__ void Min_Reduction(float *min, int *id, int offset, int gap, int NC);
__global__ void Min_Further(float *min, int *id, int NC, int BBPG, int BPG);
__global__ void Update_ID(int *id, int *id_old, int NC, int* flag);
__device__ float Calculate_Distance(float *x, int i, int j, int Dim);


void Allocate_Memory(){   
    h_flag=(int*)malloc(BPGC*sizeof(int));
    h_label=(int*)malloc(2*NP2*sizeof(int));
    h_id=(int*)malloc(NC*sizeof(int));
    h_id_old=(int*)malloc(NC*sizeof(int));
    h_x=(float*)malloc(NP*Dim*sizeof(float));
    
     
    size_t size;

    size=2*NP2*sizeof(int);
    cudaMalloc((void**) &d_label, size);
    size=BPG*NC*sizeof(int);
    cudaMalloc((void**) &d_id, size);
    size=NC*sizeof(int);
    cudaMalloc((void**) &d_id_old, size);
    size=BPG*NC*sizeof(int);
    cudaMalloc((void**) &d_nlist, size);
    size=(NC2+1)*sizeof(int);
    cudaMalloc((void**) &d_prefix, size);
    size=BPGC*sizeof(int);
    cudaMalloc((void**) &d_flag, size);
    size=NP*Dim*sizeof(float);
    cudaMalloc((void**) &d_x, size);
    size=NP*sizeof(float);
    cudaMalloc((void**) &d_sum, size);
    size=BPG*NC*sizeof(float);
    cudaMalloc((void**) &d_min, size);
    

/*
    cudaError_t Error;

    size=2*NP2*sizeof(int);
    Error=cudaMalloc((void**) &d_label, size);
    printf("CUDA Error (malloc d_label with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=BPG*NC*sizeof(int);
    Error=cudaMalloc((void**) &d_id, size);
    printf("CUDA Error (malloc d_id with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=NC*sizeof(int);
    Error=cudaMalloc((void**) &d_id_old, size);
    printf("CUDA Error (malloc d_id with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=BPG*NC*sizeof(int);
    Error=cudaMalloc((void**) &d_nlist, size);
    printf("CUDA Error (malloc d_nlist with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=(NC2+1)*sizeof(int);
    Error=cudaMalloc((void**) &d_prefix, size);
    printf("CUDA Error (malloc d_prefix with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=BPGC*sizeof(int);
    Error=cudaMalloc((void**) &d_flag, size);
    printf("CUDA Error (malloc d_flag with size %ld) = %s\n", size, cudaGetErrorString(Error));

    size=NP*Dim*sizeof(float);
    Error=cudaMalloc((void**) &d_x, size);
    printf("CUDA Error (malloc d_x with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=NP*sizeof(float);
    Error=cudaMalloc((void**) &d_sum, size);
    printf("CUDA Error (malloc d_sum with size %ld) = %s\n", size, cudaGetErrorString(Error));
    size=BPG*NC*sizeof(float);
    Error=cudaMalloc((void**) &d_min, size);
    printf("CUDA Error (malloc d_min with size %ld) = %s\n", size, cudaGetErrorString(Error));
*/
  

}
void Send_To_Device(){
    size_t size;
    size=NP*Dim*sizeof(float);
    cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);
    size=NC*sizeof(int);
    cudaMemcpy(d_id, h_id, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_id_old, h_id_old, size, cudaMemcpyHostToDevice);
/*
    cudaError_t Error;
    size=NP*Dim*sizeof(float);
    Error=cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_x -> d_x) with size %ld = %s\n",size, cudaGetErrorString(Error));
    size=NC*sizeof(int);
    Error=cudaMemcpy(d_id, h_id, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_id -> d_id) with size %ld = %s\n",size, cudaGetErrorString(Error));
    Error=cudaMemcpy(d_id_old, h_id_old, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_id_old -> d_id_old) with size %ld = %s\n",size, cudaGetErrorString(Error));
*/
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
    size_t size;

    size=NC*sizeof(int);
    cudaMemcpy(h_id, d_id, size, cudaMemcpyDeviceToHost);
    size=NP2*sizeof(int);
    cudaMemcpy(h_label, d_label, size, cudaMemcpyDeviceToHost);
    

/*
    cudaError_t Error;
    size=NC*sizeof(int);
    Error=cudaMemcpy(h_id, d_id, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_id -> h_id) with size %ld = %s\n",size, cudaGetErrorString(Error));
    size=NP2*sizeof(int);
    Error=cudaMemcpy(h_label, d_label, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_id -> h_id) with size %ld = %s\n",size, cudaGetErrorString(Error));
*/
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

/*  
    cudaError_t Error;
    Error=cudaFree(d_label);
    printf("CUDA Error ( free d_label )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_id);
    printf("CUDA Error ( free d_id )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_id_old);
    printf("CUDA Error ( free d_id_old )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_nlist);
    printf("CUDA Error ( free d_nlist )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_prefix);
    printf("CUDA Error ( free d_prefix )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_flag);
    printf("CUDA Error ( free d_flag )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_x);
    printf("CUDA Error ( free d_x )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_min);
    printf("CUDA Error ( free d_min )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_sum);
    printf("CUDA Error ( free d_sum )=%s\n", cudaGetErrorString(Error));
*/
}

char Judge(){
    
    char flag=0;
    int i;
    size_t size=sizeof(int);
    Update_ID<<<BPGC,TPB>>>(d_id, d_id_old, NC, d_flag);
    cudaMemcpy(h_flag, d_flag, size, cudaMemcpyDeviceToHost);
/*
    cudaError_t Error;
    Error=cudaMemcpy(h_flag, d_flag, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_id -> h_id) with size %ld = %s\n",size, cudaGetErrorString(Error));
*/
    for(i=0; i<BPGC; i++){
        flag|=h_flag[i];
    }
    return flag;
}


void Prefix_sum(int *d_a, int *d_b){
    printf("TPB=%d, BPGC=%d, NC2=%d, N_total=%d\n",TPB, BPGC, NC2, NP);
    Prefix_Up_Sweep<<<BPGC,TPB>>>(d_a, NC2, NP);
    Prefix_Down_Sweep<<<BPGC,TPB>>>(d_a, d_b);
    if(BPG>1){
        Offest_Between_Block<<<1,1>>>(d_b, BPG);
        Add_Offset<<<BPGC, TPB>>>(d_a, d_b);
    }
}

void Bitonic_Sort(int* d_a, char flag){
    int j,k;

//    printf("For bitonic sort, TPB=%d, BPG=%d\n",TPB, BPG2);
    
    for(k=2; k<=NP2; k<<=1){
        for(j=k>>1; j>0; j>>=1){
            Bitonic_Sort_Step<<<BPG2,TPB>>>(d_a, j, k, NP2, flag);
        }
    }
}

void Label_Point(){
    int BBPG=(int)((BPG+TPB-1)/TPB);
    int i, offset;
    
    Nearest_Medroid<<<BPG,TPB>>>(d_x, d_min, d_label,d_id,d_nlist,NC,NP, NP2,Dim);    
    //reduction
    for(i=0; i<NC; i++){
        offset=i*BPG;
        Nlist_Reduction<<<BBPG,TPB>>>(d_nlist, offset, BPG, NC);
    }
    //gather to prefix
    Nlist_To_Prefix<<<BPGC,TPB>>>(d_nlist, d_prefix, NC, NC2, BBPG, BPG, NP);
    //prefix_sum
    Prefix_sum(d_prefix, d_nlist);
    //sort
 
    Bitonic_Sort(d_label, 0);
}

void Find_Medroid(){
    int BBPG=(int)((BPG+TPB-1)/TPB);
    int i, offset; 

    Distance_Between_Cluster<<<BPG, TPB>>>(d_x, d_min, d_sum, d_label, d_prefix, d_id, Dim, NC, NP, NP2);
    for(i=0; i<NC; i++){
        offset=i*BPG;
        Min_Reduction<<<BBPG,TPB>>>(d_min, d_id, offset, BPG, NC);
    }
    
    Min_Further<<<BPGC,TPB>>>(d_min, d_id, NC, BBPG, BPG);
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

__global__ void Nearest_Medroid(float *x,float *min,int *label,int *id, int *nlist, int NC, int N, int N2, int Dim){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    float min2,temp;
    __syncthreads();
    if(i<NC*gridDim.x){
        nlist[i]=0;
        min[i]=1e15;
    }
    if(i<N){
        min2=SUM_MAX;
        //initial label for array
        label[i]=-1;
        label[i+N2]=i;
        for(int i1=0; i1<NC; i1++){
            temp=Calculate_Distance(x,i,id[i1],Dim);
            if(temp<min2){
                label[i]=i1;
                min2=temp;
            }
        }
        atomicAdd(&nlist[blockIdx.x+label[i]*gridDim.x],1);
    }
    if(i<(N2-N)){
        //initial label for redundent
        label[i+N]=NC+2;
        label[i+N+N2]=N+1;
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

__global__ void Distance_Between_Cluster(float *x, float *min, float *sum, int *label, int *prefix, int *id, int Dim, int NC, int N, int N2){
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int I = threadIdx.x;
    int up,down;
    __shared__ float temp[1024];
    __shared__ int tag[1024];
    if(i<NC*gridDim.x){
        id[i]=N+2;
    }
    if(i<N){
        sum[i]=0;
        for(int i1=prefix[label[i]]; i1<prefix[label[i]+1]; i1++ ){
            sum[i]+=Calculate_Distance(x,label[i+N2],label[i1+N2],Dim)/(prefix[label[i]+1]-prefix[label[i]]);
        }
 
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
    //no atomic min function for float, therefore doing reduction inside a block for each cluster; 

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

__global__ void Min_Reduction(float *min, int *id, int offset, int gap, int NC){
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

__global__ void Min_Further(float *min, int *id, int NC, int BBPG,int BPG){
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

__device__ float Calculate_Distance(float *x, int i, int j, int Dim){
    float temp=0;
    int k=0;
    for(k=0; k<Dim; k++){
        temp+=(x[k+i*Dim]-x[k+j*Dim])*(x[k+i*Dim]-x[k+j*Dim]);
	}
	return sqrtf(temp);
}