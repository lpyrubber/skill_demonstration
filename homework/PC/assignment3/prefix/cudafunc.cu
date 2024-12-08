#include "main.h"

void Prefix_sum();
__global__ void Prefix_Up_Sweep(int *d_a, int N, int N_total);
__global__ void Prefix_Down_Sweep(int *d_a, int *d_b);
__global__ void Offest_Between_Block(int *d_b, int N);
__global__ void Add_Offset(int *d_a, int *d_b);

void Allocate_Memory(){
    cudaError_t Error;
    size_t size=(N2+1)*sizeof(int);
    h_a=(int*)malloc(size);
    Error=cudaMalloc((void**) &d_a, size);
    printf("CUDA Error (malloc d_a) = %s\n", cudaGetErrorString(Error));
    size=(BPG+1)*sizeof(int);
    Error=cudaMalloc((void**) &d_b, size);
    printf("CUDA Error (malloc d_b) = %s\n", cudaGetErrorString(Error));

}
void Send_To_Device(){
    cudaError_t Error;
    size_t size=(N2+1)*sizeof(int);
    Error=cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_a -> d_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void GPU_Compute(){
    Prefix_sum();
}

void Send_To_Host(){
    cudaError_t Error;
    size_t size=(N2+1)*sizeof(int);
    Error=cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_a -> h_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void Free_Memory(){
    cudaError_t Error;
    free(h_a);
    Error=cudaFree(d_a);
    printf("CUDA Error ( free d_a )=%s\n", cudaGetErrorString(Error));
    Error=cudaFree(d_b);
    printf("CUDA Error ( free d_b )=%s\n", cudaGetErrorString(Error));

}

void Prefix_sum(){
    printf("TPB=%d, BPG=%d, N2=%d, N_total=%d\n",TPB, BPG, N2, N_total);
    Prefix_Up_Sweep<<<BPG,TPB>>>(d_a, N2, N_total);
    Prefix_Down_Sweep<<<BPG,TPB>>>(d_a, d_b);
    if(BPG>1){
        Offest_Between_Block<<<1,1>>>(d_b, BPG);
        Add_Offset<<<BPG, TPB>>>(d_a, d_b);
    }
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