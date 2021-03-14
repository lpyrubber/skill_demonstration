#include "main.h"

#define DEBUG 1

__global__ void GPU_scalar( float *a, float *b );

void Allocate_Memory( float **h_a , float **d_a, float **d_b ){
	cudaError_t Error;
	size_t size;
	size = N * sizeof( float );
	*h_a = ( float * )malloc( size );
	Error = cudaMalloc( ( void** )d_a, size );
	if(DEBUG) printf("CUDA Error ( malloc d_a ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_b, 2 * sizeof( float ) );
	if(DEBUG) printf("CUDA Error ( malloc d_b ) = %s\n", cudaGetErrorString( Error ) );
}

void Compute( float *a, float *b ){
	GPU_scalar<<< BPG , TPB >>>( a , b );
}

void Send_To_Host( float *h_a, float *d_a){
	cudaError_t Error;
	size_t size;
	size = N * sizeof( float );
	Error = cudaMemcpy( h_a, d_a, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error ( h_a -> d_a ) = %s\n", cudaGetErrorString( Error ) );
}

void Free_Memory( float **h_a, float **d_a, float **d_b ){
	cudaError_t Error;
	free( *h_a );
	Error = cudaFree( *d_a );
	if(DEBUG) printf( "CUDA Error ( free d_a ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_b );
	if(DEBUG) printf( "CUDA Error ( free d_b ) = %s\n", cudaGetErrorString( Error ) );
}

__global__ void GPU_scalar( float *a, float *b){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if( i < N ){
		if( i == 0 ){
			b[0] = 2.5;
		}
		__syncthreads();
		a[ i ] = b[0];
		__syncthreads();
		if( i == 0 ){
			b[ 1 ] = 2 * b[ 0 ];
		}
		if( i== 0 ){
			a[ i ] = b[ 1 ];
		}
	}
}

