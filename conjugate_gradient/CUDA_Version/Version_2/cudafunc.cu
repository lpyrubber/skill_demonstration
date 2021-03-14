#include "main.h"

__global__ void GPU_CG_1( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h);
__global__ void GPU_CG_2( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h);
__global__ void GPU_Init( float *a, float *b, float *c, float *d, float *e,  float *f, float *g );

void Allocate_Memory( float **h_a, float **h_b, float **h_c, float **d_a, float **d_b, float **d_c \
		    , float **d_d, float **d_e, float **d_f, float **d_g, float **d_h){
	size_t size;
	cudaError_t Error;

	size = N * N * sizeof(float);
	*h_a = (float*)malloc( size );
	Error = cudaMalloc( ( void** )d_a, size);
	if(DEBUG) printf("CUDA Error ( malloc d_a ) = %s\n", cudaGetErrorString( Error ) );
	size = N * sizeof(float);
	*h_b = (float*)malloc( size );
	*h_c = (float*)malloc( size );
	Error = cudaMalloc( ( void** )d_b, size);
	if(DEBUG) printf("CUDA Error ( malloc d_b ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_c, size);
	if(DEBUG) printf("CUDA Error ( malloc d_c ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_d, size);
	if(DEBUG) printf("CUDA Error ( malloc d_d ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_e, size);
	if(DEBUG) printf("CUDA Error ( malloc d_e ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_f, size);
	if(DEBUG) printf("CUDA Error ( malloc d_f ) = %s\n", cudaGetErrorString( Error ) );
	size = BPG * sizeof(float);
	Error = cudaMalloc( ( void** )d_g, size);
	if(DEBUG) printf("CUDA Error ( malloc d_g ) = %s\n", cudaGetErrorString( Error ) );
	size = 8 * sizeof(float);
	Error = cudaMalloc( ( void** )d_h, size);
	if(DEBUG) printf("CUDA Error ( malloc d_h ) = %s\n", cudaGetErrorString( Error ) );
}

void Initial( float *a, float *b, float *c ){
	int i , j;
	for( j = 0 ; j < N ; ++j ){
		c[ j ] = 0;
		b[ j ] = ( j == 0 || j == N - 1 ) ? 0 : Q/K*dx*dx ;
		for( i = 0 ; i < N ; ++i ){
			if(i == j){
				a[ i + N * j ] = 2;
			}else if( j == i - 1 && j>0 && j<N-1 ){
				a[ i + N * j ] = -1;
			}else if( j == i + 1 && j>0 && j<N-1 ){
				a[ i + N * j ] = -1;
			}else{
				a[ i + N * j ] = 0;
			}
		}
	}
}

void Send_To_Device( float *h_a, float *h_b, float *h_c, float *d_a, float *d_b , float *d_c){
	cudaError_t Error;
	size_t size;
	size = N * N * sizeof( float );
	Error = cudaMemcpy( d_a, h_a, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf("CUDA Error ( h_a -> d_a ) = %s\n", cudaGetErrorString( Error ) );
	size = N * sizeof( float );
	Error = cudaMemcpy( d_b, h_b, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf("CUDA Error ( h_b -> d_b ) = %s\n", cudaGetErrorString( Error ) );
	Error = cudaMemcpy( d_c, h_c, size, cudaMemcpyHostToDevice );
	if(DEBUG) printf("CUDA Error ( h_c -> d_c ) = %s\n", cudaGetErrorString( Error ) );

}

void Device_Initial( float *a, float *b, float *c, float *d, float *e, float *f, float *g){
	printf("BPG = %d, TPB = %d\n", BPG, TPB);
	GPU_Init<<< BPG, TPB >>>( a, b, c, d, e, f, g );
}	

void Conjugate_Gradient( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h, float *temp){
	cudaError_t Error;
	size_t size;
//	printf("BPG = %d, TPB = %d\n", BPG, TPB);
	GPU_CG_1<<< BPG, TPB >>>( a, b, c, d, e, f, g, h );
	GPU_CG_2<<< BPG, TPB >>>( a, b, c, d, e, f, g, h );
	size = sizeof(float);
	Error = cudaMemcpy( temp, &h[0], size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf("CUDA Error ( h_err -> d_err ) = %s\n", cudaGetErrorString( Error ) );
}

void Send_To_Host( float *h_a, float *d_a ){
	cudaError_t Error;
	size_t size;
	size = N * sizeof( float );
	Error = cudaMemcpy( h_a, d_a, size, cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error ( d_a -> h_a ) = %s\n", cudaGetErrorString( Error ) );
}

void Save_Result(float *a){
	FILE *fp;
	int i;
	fp = fopen( "data.txt", "w");
	for( i = 0 ; i < N ; ++i ){
		fprintf(fp, "%f\n", a[i]);
	}
	fclose(fp);
}

void Free_Memory( float **h_a, float **h_b, float **h_c, float **d_a, float **d_b, float **d_c \
		, float **d_d, float **d_e, float **d_f, float **d_g, float **d_h){
	cudaError_t Error;
	free( *h_a );
	free( *h_b );
	free( *h_c );
	Error = cudaFree( *d_a );
	if(DEBUG) printf("CUDA Error ( free d_a )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_b );
	if(DEBUG) printf("CUDA Error ( free d_b )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_c );
	if(DEBUG) printf("CUDA Error ( free d_c )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_d );
	if(DEBUG) printf("CUDA Error ( free d_d )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_e );
	if(DEBUG) printf("CUDA Error ( free d_e )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_f );
	if(DEBUG) printf("CUDA Error ( free d_f )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_g );
	if(DEBUG) printf("CUDA Error ( free d_g )=%s\n", cudaGetErrorString(Error));
	Error = cudaFree( *d_h );
	if(DEBUG) printf("CUDA Error ( free d_h )=%s\n", cudaGetErrorString(Error));
}

__global__ void GPU_Init( float *a, float *b, float *c, float *d, float *e, float *f, float *g ){
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int I = threadIdx.x;
	__shared__ float temp[TPB];
	if( i < N ){
		d[ i ] = b[ i ]; //r[ i ] = b[ i ] 
		for( int i1 = 0 ; i1 < N ; ++i1 ){
			d[ i ] -= a[ i1 + N * i ] * c[ i1 ]; //r[ i ] = A[ i ][ j ] * x[ j ]
		}
		__syncthreads();
		e[ i ] = d[ i ]; //p[ i ] = r[ i ]

		temp[ I ] = 0;
		temp[ I ] = d[ i ]*d[ i ];
		__syncthreads();
		for(int stride = blockDim.x/2; stride > 0; stride = stride/2){
			if(I<stride){
				temp[ I ] += temp[ I + stride ];
			}
			__syncthreads();
		}
		if( I == 0 ) f[ blockIdx.x ] = temp[ 0 ];
		__syncthreads();
		if( i == 0 ){
			g[ 0 ] = 0;
			for( int i1 = 0 ; i1 < gridDim.x ; ++i1){
				g[ 0 ] += f[ i1 ];
			}
		}
		__syncthreads();
	}
}

__global__ void GPU_CG_1( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h ){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int I = threadIdx.x;
	__shared__ float temp[TPB];
	if( i < N ){
		e[ i ] = 0; //s[ i ] = 0;
		for( int i1 = 0 ; i1 < N ; ++i1 ){
                        e[ i ] += a[ i1 + N * i ] * d[ i1 ]; //s[ i ] = A[ i ][ j ] * p[ j ]
                }
                __syncthreads();

                temp[ I ] = 0;
                temp[ I ] = e[ i ]*d[ i ]; //temp = s * p'
                __syncthreads();
                for(int stride = blockDim.x/2; stride > 0; stride = stride/2){
                        if(I<stride){
                                temp[ I ] += temp[ I + stride ];
                        }
                        __syncthreads();
                }
                if( I == 0 ) g[ blockIdx.x ] = temp[ 0 ];
                __syncthreads();
                if( i == 0 ){
			h[ 2 ] = 0;
                        for( int i1 = 0 ; i1 < gridDim.x ; i1++){
                                h[ 2 ] += g[ i1 ];
                        }
			h[ 3 ] = h[ 0 ] / h[ 2 ]; //alpha = rsold / temp;
		}
		__syncthreads();
	}
}
__global__ void GPU_CG_2( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h ){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int I = threadIdx.x;
	__shared__ float temp[TPB];
	if( i < N ){
		c[ i ] += h[3] * d[ i ]; //x[ i ] += alpha * p[ i ]
		f[ i ] -= h[3] * e[ i ]; //r[ i ] -= alpha * s[ i ]
		__syncthreads();
                temp[ I ] = 0;
                temp[ I ] = f[ i ]*f[ i ]; //rsnew = r * r'
                __syncthreads();
                for(int stride = blockDim.x/2; stride > 0; stride = stride/2){
                        if(I<stride){
                                temp[ I ] += temp[ I + stride ];
                        }
                        __syncthreads();
                }
                if( I == 0 ) g[ blockIdx.x ] = temp[ 0 ];
                __syncthreads();
                if( i == 0 ){
			h[ 1 ] = 0;
                        for( int i1 = 0 ; i1 < gridDim.x ; i1 ++){
                                h[ 1 ] += g[ i1 ];
                        }
               }
		__syncthreads();
		d[ i ] = f[ i ] + ( h[1] / h[0] ) * d[ i ]; // p[ i ] = r [ i ] + (rsnew/rsold) * p[ i ]
		__syncthreads();
		if( i == 0 ) h[ 0 ] = h[ 1 ];
	}
}
