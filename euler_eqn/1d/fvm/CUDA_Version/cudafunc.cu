#include "main.h"

__global__ void GPU_Calc( float *a, float *b, float *c, float *d, float *e, float *f );

void Allocate_Memory( float **h_a, float **h_b, float **h_c, float **h_d\
		   , float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f ){
	cudaError_t Error;
	size_t size;
	printf( "Allocating Memory....\n" );
	size = N * sizeof( float );
	*h_a = ( float* )malloc( 3 * size );
	*h_b = ( float* )malloc( size );
	*h_c = ( float* )malloc( size );
	*h_d = ( float* )malloc( size );

	Error = cudaMalloc( ( void** )d_a , 3 * size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_a) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_b , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_b) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_c , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_c) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_d , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_d) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_e , 3 * size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_e) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_f , 3 * size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_f) = %s\n", cudaGetErrorString( Error ) );
}

void Send_To_Device( float *h_a, float *h_b, float *h_c, float *h_d\
		   , float *d_a, float *d_b, float *d_c, float *d_d ){
	cudaError_t Error;
	size_t size;
	printf( "Sending data to Device(GPU)...\n" );
	size = N * sizeof( float );
	Error = cudaMemcpy( d_a, h_a, 3 * size, cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA ERROR Memcpy h_a -> d_a = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_b, h_b, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA ERROR Memcpy h_b -> d_b = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_c, h_c, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA ERROR Memcpy h_c -> d_c = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_d, h_d, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA ERROR Memcpy h_d -> d_d = %s\n", cudaGetErrorString( Error ) );
}

void GPU_Compute( float *d_a, float *d_b, float *d_c, float *d_d, float *d_e, float *d_f ){
	GPU_Calc<<< BPG, TPB >>>( d_a, d_b, d_c, d_d, d_e, d_f );
}

__global__ void GPU_Calc( float *a, float *b, float *c, float *d, float *e, float *f ){
	int i;
	float acc, F1, F2, F3, Fr;
	float FL1, FL2, FL3, FR1, FR2, FR3;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	if( i < N ){
		acc = sqrt( gamma * R * d[ i ] );
		Fr = c[ i ] / acc;
		if ( Fr > 1  ) Fr =  1;
		if ( Fr < -1 ) Fr = -1;
		F1 = a[ i ] * c[ i ];
		F2 = a[ i ] * c[ i ] * c[ i ] +  b[ i ] * R * d[ i ];
		F3 = c[ i ] * ( a[ i + N * 2 ] + b[ i ] * R * d[ i ] );
		f[ i         ] =  0.5*(F1*(Fr+1)+a[ i         ]*acc*(1-Fr*Fr));
		f[ i + N     ] =  0.5*(F2*(Fr+1)+a[ i + N     ]*acc*(1-Fr*Fr));
		f[ i + N * 2 ] =  0.5*(F3*(Fr+1)+a[ i + N * 2 ]*acc*(1-Fr*Fr));
		e[ i         ] = -0.5*(F1*(Fr-1)+a[ i         ]*acc*(1-Fr*Fr));
		e[ i + N     ] = -0.5*(F2*(Fr-1)+a[ i + N     ]*acc*(1-Fr*Fr));
		e[ i + N * 2 ] = -0.5*(F3*(Fr-1)+a[ i + N * 2 ]*acc*(1-Fr*Fr));
	}
	__syncthreads();
	if( i > 0 && i < N - 1 ){
		FL1 = f[ i - 1         ]+e[ i             ];
		FR1 = f[ i             ]+e[ i + 1         ];
		FL2 = f[ i - 1 + N     ]+e[ i     + N     ];
		FR2 = f[ i     + N     ]+e[ i + 1 + N     ];
		FL3 = f[ i - 1 + N * 2 ]+e[ i     + N * 2 ];
		FR3 = f[ i     + N * 2 ]+e[ i + 1 + N * 2 ];
		a[ i         ] = a[ i         ] - Z*(FR1-FL1);
		a[ i + N     ] = a[ i + N     ] - Z*(FR2-FL2);
		a[ i + N * 2 ] = a[ i + N * 2 ] - Z*(FR3-FL3);
	}
	__syncthreads();
	if( i == 0 ){
		a[ i         ] =  a[ i + 1         ];
		a[ i + N     ] = -a[ i + 1 + N     ];
		a[ i + N * 2 ] =  a[ i + 1 + N * 2 ];
	}
	if( i == N - 1 ){	
		a[ i         ] =  a[ i - 1         ];
		a[ i + N     ] = -a[ i - 1 + N     ];
		a[ i + N * 2 ] =  a[ i - 1 + N * 2 ];
	}
	__syncthreads();
	if( i < N ){	
		b[ i ] = a[ i ];
		c[ i ] = a[ i + N ]/a[i];
		d[ i ] = ( ( a[ i + N * 2 ]/a[ i ])-0.5*c[ i ]*c[ i ])/CV;
	}
}

void Send_To_Host( float *h_a, float *h_b, float *h_c, float *h_d\
		 , float *d_a, float *d_b, float *d_c, float *d_d ){
	cudaError_t Error;
	size_t size;
	size = N * sizeof( float );
	printf( "Sending data to host(GPU)...\n" );
	Error = cudaMemcpy( h_a, d_a, 3 * size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_a -> h_a = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( h_b, d_b, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_b -> h_b = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( h_c, d_c, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_c -> h_c = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( h_d, d_d, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_d -> h_d = %s\n", cudaGetErrorString( Error ) );
}

void Free( float **h_a, float **h_b, float **h_c, float **h_d\
	 , float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f ){
	cudaError_t Error;
	printf( "Free memory...\n" );
	free( *h_a );
	free( *h_b );
	free( *h_c );
	free( *h_d );
	Error = cudaFree( *d_a );
	if(DEBUG) printf( "CUDA Error cudaFree( d_a ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_b );
	if(DEBUG) printf( "CUDA Error cudaFree( d_b ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_c );
	if(DEBUG) printf( "CUDA Error cudaFree( d_c ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_d );
	if(DEBUG) printf( "CUDA Error cudaFree( d_d ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_e );
	if(DEBUG) printf( "CUDA Error cudaFree( d_e ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_f );
	if(DEBUG) printf( "CUDA Error cudaFree( d_f ) = %s\n", cudaGetErrorString( Error ) );
}
