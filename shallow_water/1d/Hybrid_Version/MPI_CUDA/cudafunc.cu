#include "main.h"

__global__ void GPU_Calc( float *d_a , float *d_b , float *d_c , int ele );

void Allocate_Memory( int id , int max , float **h_a , float **h_b , float **d_a , float **d_b , float **d_c ){
	//define variable in function
	int averow , extra , rows;
	size_t size;
	cudaError_t Error;
	if( id == MASTER ){
		*h_a = ( float* )malloc( N * sizeof( float ) );
		*h_b = ( float* )malloc( N * sizeof( float ) );
	}else{
		averow = N / ( max - 1 );
		extra  = N % ( max - 1 );
		rows = ( id <= extra ) ?  averow + 3 : averow + 2 ; // +2 for ghost blocks
		size = rows * sizeof( float ) ;
		*h_a = ( float* )malloc( size );
		*h_b = ( float* )malloc( size );
		Error = cudaSetDevice( id - 1 );
		if(DEBUG) printf( "CUDA Error cudaSetDevice( %d ) = %s\n", id - 1 , cudaGetErrorString( Error ) );
		Error = cudaMalloc( ( void** )d_a , 1 * size );
		if(DEBUG) printf( "CUDA Error cudaMalloc(d_a) = %s\n" , cudaGetErrorString( Error ) );
		Error = cudaMalloc( ( void** )d_b , 1 * size );
		if(DEBUG) printf( "CUDA Error cudaMalloc(d_b) = %s\n" , cudaGetErrorString( Error ) );
		Error = cudaMalloc( ( void** )d_c , 4 * size );
		if(DEBUG) printf( "CUDA Error cudaMalloc(d_c) = %s\n" , cudaGetErrorString( Error ) );
	}
}

void Sent_To_Device( int ele , float *h_a , float *h_b , float *d_a , float *d_b ){
	cudaError_t Error;
	size_t size;
//	printf( "Sending data to Device(GPU)....\n" );
	size = ( ele + 2 ) * sizeof( float );
	Error = cudaMemcpy( d_a , h_a , size , cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy h_a->d_a = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_b , h_b , size , cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy h_b->d_a = %s\n" , cudaGetErrorString( Error ) );
}

void Compute( int ele , float *d_a , float *d_b , float *d_c ){
	GPU_Calc<<< BPG , TPB >>>( d_a , d_b , d_c , ele);
}

__global__ void GPU_Calc( float *d_a , float *d_b , float *d_c , int ele ){
	int i;
	float vel , acc , F1  , F2  , Fr \
	    , FL1 , FL2 , FR1 , FR2 ;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	//compute flux
	if(i < ele + 2){
		vel = d_b[ i ] / d_a[ i ];
		acc = sqrt( G * d_a[ i ] );
		Fr = vel / acc ;
		if ( Fr >  1 ) Fr =  1 ;
		if ( Fr < -1 ) Fr = -1 ;
		F1 = d_a[ i ] * vel ;
		F2 = d_a[ i ] * vel * vel + 0.5 * G * d_a[ i ] * d_a[ i ] ;
		d_c[ i                 ] =   0.5 * ( F1 * ( Fr + 1 ) + d_a[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fp_a
		d_c[ i + ( ele + 2 )     ] =   0.5 * ( F2 * ( Fr + 1 ) + d_b[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fp_b
		d_c[ i + ( ele + 2 ) * 2 ] = - 0.5 * ( F1 * ( Fr - 1 ) + d_a[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fm_a
		d_c[ i + ( ele + 2 ) * 3 ] = - 0.5 * ( F2 * ( Fr - 1 ) + d_b[ i ] * acc * ( 1 - Fr * Fr ) ) ; //fm_b
	}
	__syncthreads();
	//computee u
	if( ( i > 0 ) && ( i <= ele ) ){
		FL1 = d_c[ i - 1             ] + d_c[ i     + ( ele + 2 ) * 2 ];
		FR1 = d_c[ i                 ] + d_c[ i + 1 + ( ele + 2 ) * 2 ];
		FL2 = d_c[ i - 1 + ( ele + 2 ) ] + d_c[ i     + ( ele + 2 ) * 3 ];
		FR2 = d_c[ i     + ( ele + 2 ) ] + d_c[ i + 1 + ( ele + 2 ) * 3 ];
		d_a[ i ] = d_a[ i ] - Z * ( FR1 - FL1 ) ;
		d_b[ i ] = d_b[ i ] - Z * ( FR2 - FL2 ) ;
	}
}

void Sent_To_Host( int ele , float *h_a , float *h_b , float *d_a , float *d_b ){
	cudaError_t Error;
	size_t size;
//	printf( "Sending data to Host(CPU)....\n" );
	size = ( ele + 2 ) * sizeof( float );
	Error = cudaMemcpy( h_a , d_a , size , cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_a->h_a = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMemcpy( h_b , d_b , size , cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_b->h_b = %s\n" , cudaGetErrorString( Error ) );
}

void Free_Memory( int id , float **h_a , float **h_b , float **d_a , float **d_b , float **d_c ){
	cudaError_t Error;
	printf( "Free memory....\n" );
	if(id == MASTER){
		free( *h_a );
		free( *h_b );
	}else{
		free( *h_a );
		free( *h_b );
		Error = cudaFree( *d_a );
		if(DEBUG) printf( "CUDA Error cudaFree(d_a) = %s\n" , cudaGetErrorString( Error ) );
		Error = cudaFree( *d_b );
		if(DEBUG) printf( "CUDA Error cudaFree(d_b) = %s\n" , cudaGetErrorString( Error ) );
		Error = cudaFree( *d_c );
		if(DEBUG) printf( "CUDA Error cudaFree(d_c) = %s\n" , cudaGetErrorString( Error ) );
	}
}
