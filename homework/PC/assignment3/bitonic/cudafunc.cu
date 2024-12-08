#include "main.h"

void Bitonic_Sort();
__global__ void Bitonic_Sort_Step(int *d_a, int j, int k, int N, char flag);


void Allocate_Memory(){
    cudaError_t Error;
    size_t size=2*N2*sizeof(int);
    h_a=(int*)malloc(size);
    printf("before, d_a address=%p\n",d_a);
    Error=cudaMalloc((void**) &d_a, size);
    printf("CUDA Error (malloc d_a) = %s\n", cudaGetErrorString(Error));
    printf("after d_a address= %p\n",d_a);
}
void Send_To_Device(){
    cudaError_t Error;
    size_t size=2*N2*sizeof(int);
    Error=cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
    printf("CUDA Error (h_a -> d_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void GPU_Compute(){
    Bitonic_Sort();
}
void Send_To_Host(){
    cudaError_t Error;
    size_t size=2*N2*sizeof(int);
    printf("h_a = %p\n", h_a);
    Error=cudaMemcpy(h_a, d_a, size, cudaMemcpyDeviceToHost);
    printf("CUDA Error (d_a -> h_a) with size %ld = %s\n",size, cudaGetErrorString(Error));
}
void Free_Memory(){
    cudaError_t Error;
    free(h_a);
    Error=cudaFree(d_a);
    printf("CUDA Error ( free d_a )=%s\n", cudaGetErrorString(Error));
}

void Bitonic_Sort(){
    int j,k;
    k=2;
    j=1;
    printf("TPB=%d, BPG=%d, j=%d, k=%d d_a=%p\n",TPB, BPG, j, k, d_a);
    
    for(k=2; k<=N2; k<<=1){
        for(j=k>>1; j>0; j>>=1){
            Bitonic_Sort_Step<<<BPG,TPB>>>(d_a, j, k, N2, 1);
        }
    }

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