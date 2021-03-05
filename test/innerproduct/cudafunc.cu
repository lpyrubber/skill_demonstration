#include "main.h"

__global__ void GPU_reduct( float *a, float *b);

void Allocate_Memory(float **h_a, float **h_b, float **d_a, float **d_b){
	cudaError_t Error;
	size_t size;
	size = N * sizeof(float);
	*h_a = (float*)malloc(size);
	Error = cudaMalloc( ( void** )d_a, size);
	if(DEBUG) printf("CUDA Error ( malloc d_a ) = %s\n", cudaGetErrorString( Error ) );
	size = BPG * sizeof(float);
	*h_b = (float*)malloc(size);
	Error = cudaMalloc( ( void** )d_b, size);
	if(DEBUG) printf("CUDA Error ( malloc d_b ) = %s\n", cudaGetErrorString( Error ) );
}

void Initial(float *h_a ){
	int i;
	for(i = 0 ; i < N ; i++){
		h_a[ i ] = 2.5;
	}
}
void Send_To_Device( float *h_a , float *d_a ){
	cudaError_t Error;
	size_t size;
	size = N * sizeof( float );
	Error = cudaMemcpy( d_a, h_a, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf("CUDA Error ( h_a -> d_a ) = %s\n", cudaGetErrorString( Error ) );
}
void GPU_Dot( float *a, float *b ){
	GPU_reduct<<< BPG, TPB >>>(a, b);
}
void Send_To_Host( float *h_a, float *d_a ){
	cudaError_t Error;
	size_t size;
	size = BPG * sizeof( float );
	Error = cudaMemcpy( h_a, d_a, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error ( d_a -> h_a ) = %s\n", cudaGetErrorString( Error ) );
}
void Free_Memory( float **h_a, float **h_b, float **d_a, float **d_b ){
	cudaError_t Error;
	free( *h_a );
	free( *h_b );
	Error = cudaFree( *d_a );
	if(DEBUG) printf("CUDA Error ( free d_a )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_b );
	if(DEBUG) printf("CUDA Error ( free d_b )=%s\n", cudaGetErrorString(Error));
}

__global__ void GPU_reduct( float *a ,float *b){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int I = threadIdx.x;
	__shared__ float temp[TPB];
	if( i < N ){
		temp[ I ] = 0;
		temp[ I ] = a[ i ]*a[ i ];
		__syncthreads();
		for(int stride = blockDim.x/2; stride > 0; stride >>= 1){
			if(I<stride){
				temp[ I ] += temp[ I + stride ];
			}
			__syncthreads();
		}
		__syncthreads();
		if( I == 0 ) b[ blockIdx.x ] = temp[ 0 ];
		__syncthreads();
//		if( I == 0 ) atomicAdd( b , temp[ 0 ] );
		if( i == 0 ){
			float sum = 0;
			for( int i1 = 0 ; i1 < gridDim.x ; i1 += 2){
				sum += b[ i1 ];
			}
			b[ 0 ] = sum;
		}else if( i == 1 ){
			float sum = 0;
			for( int i1 = 1 ; i1 < gridDim.x ; i1 += 2){
				sum += b[ i1 ];
			}
			b[ 1 ] = sum;
		}
		__syncthreads();
		if( i == 0 ){
			b[ 0 ] += b[ 1 ];
		}
	}
}
