#include "main.h"

__global__ void GPU_Calc( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h, float *i );

void Allocate_Memory( float **h_a, float **h_b, float **h_c, float **h_d\
		   , float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i ){
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
	Error = cudaMalloc( ( void** )d_f , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_f) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_g , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_g) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_h , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_h) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_i , size );
	if(DEBUG) printf( "CUDA ERROR cudaMalloc(d_i) = %s\n", cudaGetErrorString( Error ) );
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

void GPU_Compute( float *d_a, float *d_b, float *d_c, float *d_d, float *d_e, float *d_f, float *d_g, float *d_h, float *d_i ){
	GPU_Calc<<< BPG, TPB >>>( d_a, d_b, d_c, d_d, d_e, d_f, d_g, d_h, d_i );
}

__global__ void GPU_Calc( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h, float *i ){
	int id;
	float acc, M, Mt, P;
	float FL[3], FR[3];
	id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id<N){
		acc=sqrt(gamma*R*d[id]);
		M=c[id]/acc;
		P=b[id]*R*d[id];
		Mt=fabs(M);
		if(Mt<=1){
			f[id]= 0.25*P*(M+1)*(M+1)*(2-M);
			g[id]= 0.25*P*(M-1)*(M-1)*(2+M);
			h[id]= 0.25*(M+1)*(M+1);
			i[id]=-0.25*(M-1)*(M-1);
		}else{
			f[id]= 0.5*P*(M+Mt)/M;
			g[id]= 0.5*P*(M-Mt)/M;
			h[id]= 0.5*(M+Mt);
			i[id]=-0.5*(M-Mt);
		}
		e[id]=b[id]*acc;
		e[id+N]=b[id]*acc*c[id];
		e[id+N*2]=b[id]*acc*(CP*d[id]+0.5*c[id]*c[id]);
	}
	__syncthreads();
	if(id>1&&id<N-1){
		Mt=h[id-1]+i[id];
		FL[0]=Mt*(e[id-1    ]+e[id    ])*0.5-fabs(Mt)*(e[id    ]-e[id-1    ])*0.5;
		FL[1]=Mt*(e[id-1+N  ]+e[id+N  ])*0.5-fabs(Mt)*(e[id+N  ]-e[id-1+N  ])*0.5+(f[id-1]+g[id]);
		FL[2]=Mt*(e[id-1+N*2]+e[id+N*2])*0.5-fabs(Mt)*(e[id+N*2]-e[id-1+N*2])*0.5;
		Mt=h[id]+i[id+1];
		FR[0]=Mt*(e[id    ]+e[id+1    ])*0.5-fabs(Mt)*(e[id+1    ]-e[id    ])*0.5;
		FR[1]=Mt*(e[id+N  ]+e[id+1+N  ])*0.5-fabs(Mt)*(e[id+1+N  ]-e[id+N  ])*0.5+(f[id]+g[id+1]);
		FR[2]=Mt*(e[id+N*2]+e[id+1+N*2])*0.5-fabs(Mt)*(e[id+1+N*2]-e[id+N*2])*0.5;
		a[id    ]-=Z*(FR[0]-FL[0]);
		a[id+N  ]-=Z*(FR[1]-FL[1]);
		a[id+N*2]-=Z*(FR[2]-FL[2]);
	}
	__syncthreads();
	if( id == 0 ){
		a[ id         ] =  a[ id + 1         ];
		a[ id + N     ] = -a[ id + 1 + N     ];
		a[ id + N * 2 ] =  a[ id + 1 + N * 2 ];
	}
	if( id == N - 1 ){	
		a[ id         ] =  a[ id - 1         ];
		a[ id + N     ] = -a[ id - 1 + N     ];
		a[ id + N * 2 ] =  a[ id - 1 + N * 2 ];
	}
	__syncthreads();
	if( id < N ){	
		b[ id ] = a[ id ];
		c[ id ] = a[ id + N ]/a[id];
		d[ id ] = ( ( a[ id + N * 2 ]/a[ id ])-0.5*c[ id ]*c[ id ])/CV;
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
	 , float **d_a, float **d_b, float **d_c, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h, float **d_i ){
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
	Error = cudaFree( *d_g );
	if(DEBUG) printf( "CUDA Error cudaFree( d_g ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_h );
	if(DEBUG) printf( "CUDA Error cudaFree( d_h ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaFree( *d_i );
	if(DEBUG) printf( "CUDA Error cudaFree( d_i ) = %s\n", cudaGetErrorString( Error ) );
}
