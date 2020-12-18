#include "main.h"

__global__ void GPU_Calc( float *d_a , float *d_b, float *d_c );

void Allocate_Memory( float **h_a , float **h_b , float **d_a , float **d_b , float **d_c ){
	cudaError_t Error;
	size_t size;
	printf( "Allocating memory....\n" );
	size = ( N + 2 ) * sizeof( float );
	*h_a = ( float* )malloc( size );
	*h_b = ( float* )malloc( size );
	Error = cudaMalloc( ( void** )d_a , 2 * size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_a) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_b , 2 * size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_b) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_c , 4 * size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_c) = %s\n" , cudaGetErrorString( Error ) );
	
}

void Sent_To_Device( float *h_a , float *h_b , float *d_a , float *d_b ){
	cudaError_t Error;
	size_t size;
	printf( "Sending data to Device(GPU)....\n" );
	size = ( N + 2 ) * sizeof( float );
	Error = cudaMemcpy( d_a , h_a , size , cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy h_a->d_a = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_b , h_b , size , cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy h_b->d_a = %s\n" , cudaGetErrorString( Error ) );
}

void GPU_Compute( float *d_a , float *d_b , float *d_c ){
	cudaError_t Error;
	size_t size ;
	size = ( N + 2 ) * sizeof( float ) ; 
	GPU_Calc<<< BPG , TPB >>>( d_a , d_b , d_c );
	Error = cudaMemcpy( d_a , d_a + N + 2 , size , cudaMemcpyDeviceToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy d_a+N+2->d_a = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_b , d_b + N + 2 , size , cudaMemcpyDeviceToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy d_b+N+2->d_b = %s\n" , cudaGetErrorString( Error ) );	 
}

__global__ void GPU_Calc( float *d_a , float *d_b , float *d_c ){
	int i;
	float vel , acc , F1  , F2  , Fr \
	    , FL1 , FL2 , FR1 , FR2 ;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	//compute flux
	if(i < N + 2){
		vel = d_b[ i ] / d_a[ i ];
		acc = sqrt( G * d_a[ i ] );
		Fr = vel / acc ;
		if ( Fr >  1 ) Fr =  1 ;
		if ( Fr < -1 ) Fr = -1 ;
		F1 = d_a[ i ] * vel ;
		F2 = d_a[ i ] * vel * vel + 0.5 * G * d_a[ i ] * d_a[ i ] ;
		d_c[ i                 ] =   0.5 * ( F1 * ( Fr + 1 ) + d_a[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fp_a
		d_c[ i + ( N + 2 )     ] =   0.5 * ( F2 * ( Fr + 1 ) + d_b[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fp_b
		d_c[ i + ( N + 2 ) * 2 ] = - 0.5 * ( F1 * ( Fr - 1 ) + d_a[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fm_a
		d_c[ i + ( N + 2 ) * 3 ] = - 0.5 * ( F2 * ( Fr - 1 ) + d_b[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fm_b
	}
	__syncthreads();
	//computee u
	if( ( i > 0 ) && ( i <= N ) ){
		FL1 = d_c[ i - 1             ] + d_c[ i     + ( N + 2 ) * 2 ];
		FR1 = d_c[ i                 ] + d_c[ i + 1 + ( N + 2 ) * 2 ];
		FL2 = d_c[ i - 1 + ( N + 2 ) ] + d_c[ i     + ( N + 2 ) * 3 ];
		FR2 = d_c[ i     + ( N + 2 ) ] + d_c[ i + 1 + ( N + 2 ) * 3 ];
		d_a[ i + N + 2 ] = d_a[ i ] - Z * ( FR1 - FL1 ) ;
		d_b[ i + N + 2 ] = d_b[ i ] - Z * ( FR2 - FL2 ) ;
	}
	//wait for a while
	__syncthreads();
	//set reflective boundary
	if(  i == 0 ){
		d_a[ i + N + 2] =   d_a[ i + 1 ];
		d_b[ i + N + 2] = - d_b[ i + 1 ];
	}
	if( i == N + 1 ){
		d_a[ i + N + 2] =   d_a[ i - 1 ];
		d_b[ i + N + 2] = - d_b[ i - 1 ];
	}
}

void Sent_To_Host( float *h_a , float *h_b , float *d_a , float *d_b ){
	cudaError_t Error;
	size_t size;
	printf( "Sending data to Host(CPU)....\n" );
	size = ( N + 2 ) * sizeof( float );
	Error = cudaMemcpy( h_a , d_a , size , cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_a->h_a = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMemcpy( h_b , d_b , size , cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_b->h_b = %s\n" , cudaGetErrorString( Error ) );
}

void Free( float **h_a , float **h_b , float **d_a , float **d_b , float **d_c ){
	cudaError_t Error;
	printf( "Free memory....\n" );
	free( *h_a );
	free( *h_b );
	Error = cudaFree( *d_a );
	if(DEBUG) printf( "CUDA Error cudaFree(d_a) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaFree( *d_b );
	if(DEBUG) printf( "CUDA Error cudaFree(d_b) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaFree( *d_c );
	if(DEBUG) printf( "CUDA Error cudaFree(d_c) = %s\n" , cudaGetErrorString( Error ) );
}
