#include "main.h"

__global__ void GPU_Calc( float *a , float *b, float *c );

void Allocate_Memory( float **h_a , float **d_a , float **d_b , float **d_c ){
	cudaError_t Error;
	size_t size;
	printf( "Allocating memory....\n" );
	size =  3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
	*h_a = ( float* )malloc( size );
	Error = cudaMalloc( ( void** )d_a ,     size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_a) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_b , 2 * size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_b) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaMalloc( ( void** )d_c , 2 * size );
	if(DEBUG) printf( "CUDA Error cudaMalloc(d_c) = %s\n" , cudaGetErrorString( Error ) );
	
}

void Sent_To_Device( float *h_a , float *d_a ){
	cudaError_t Error;
	size_t size;
	printf( "Sending data to Device(GPU)....\n" );	
	size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
	Error = cudaMemcpy( d_a , h_a , size , cudaMemcpyHostToDevice );
	if(DEBUG) printf( "CUDA Error Memcpy h_a->d_a = %s\n" , cudaGetErrorString( Error ) );
}

void GPU_Compute( float *d_a , float *d_b , float *d_c ){
	GPU_Calc<<< BPG , TPB >>>( d_a , d_b , d_c ); 
}

__global__ void GPU_Calc( float *d_a , float *d_b , float *d_c ){
	int i;
	float velx , vely , acc , F1  , F2 ,  F3 , Frx, Fry \
	    , FL1x , FL2x , FL3x , FR1x , FR2x , FR3x \
	    , FL1y , FL2y , FL3y , FR1y , FR2y , FR3y;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	//compute flux
	if( i < ( Nx + 2 ) * ( Ny + 2 ) ){
		velx = d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] / d_a[ i ];
		vely = d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] / d_a[ i ];
		acc    = sqrt( G * d_a[ i  ] );
		Frx = velx / acc;
		Fry = vely / acc;
		if (Frx >  1 ) Frx =  1;
		if (Frx < -1 ) Frx = -1;
		if (Fry >  1 ) Fry =  1;
		if (Fry < -1 ) Fry = -1;
		//x-dir
		F1 = d_a[ i ] * velx;
		F2 = d_a[ i ] * velx * velx\
		   + 0.5 * G * d_a[ i ] * d_a[ i ];
		F3 = d_a[ i ] * velx * vely;
		d_b[ i                                ]\
			=  0.5 * ( F1 * ( Frx + 1 )\
			+ d_a[ i                               ] * acc * ( 1 - Frx * Frx ) );
		d_b[ i  + ( Nx + 2 ) * ( Ny + 2 )     ]\
			=  0.5 * ( F2 * ( Frx + 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Frx * Frx ) );
		d_b[ i  + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			=  0.5 * ( F3 * ( Frx + 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
		d_c[ i                                ]\
			= -0.5 * ( F1 * ( Frx - 1 )\
			+ d_a[ i                               ] * acc * ( 1 - Frx * Frx ) );
		d_c[ i + ( Nx + 2 ) * ( Ny + 2 )     ]\
			= -0.5 * ( F2 * ( Frx - 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Frx * Frx ) );
 		d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]
			= -0.5 * ( F3 * ( Frx - 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
		//y-dir
		F1 = d_a[ i ] * vely;
		F2 = d_a[ i ] * vely * velx;
		F3 = d_a[ i ] * vely * vely\
		   + 0.5 * G * d_a[ i ] * d_a[ i ];
		d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			=  0.5 * ( F1 * ( Fry + 1 )\
			+ d_a[ i                               ] * acc * ( 1 - Fry * Fry ) );
		d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			=  0.5 * ( F2 * ( Fry + 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			=  0.5 * ( F3 * ( Fry + 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );
		d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			= -0.5 * ( F1 * ( Fry - 1 )\
			+ d_a[ i                               ] * acc * ( 1 - Fry * Fry ) );
		d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			= -0.5 * ( F2 * ( Fry - 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			= -0.5 * ( F3 * ( Fry - 1 )\
			+ d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );	
	}
	__syncthreads();
	//computee u
	if( ( i / ( Nx + 2 ) > 0 ) && ( i / ( Nx + 2 ) <= Nx ) && ( i % ( Nx + 2 ) > 0 ) && ( i % ( Nx + 2 ) <= Nx ) ){
			FL1x = d_b[ i - 1                                ]\
			     + d_c[ i                                    ];
			FR1x = d_b[ i                                    ]\
			     + d_c[ i + 1                                ];
			FL2x = d_b[ i - 1 + ( Nx + 2 ) * ( Ny + 2 )      ]\
			     + d_c[ i + ( Nx + 2 ) * ( Ny + 2 )          ];
			FR2x = d_b[ i + ( Nx + 2 ) * ( Ny + 2 )          ]\
			     + d_c[ i + 1 + ( Nx + 2 ) * ( Ny + 2 )      ];
			FL3x = d_b[ i - 1 + ( Nx + 2 ) * ( Ny +  2 ) * 2 ]\
			     + d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2      ];
			FR3x = d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2      ]\
			     + d_c[ i + 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2  ];

			FL1y = d_b[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			     + d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3              ];
			FR1y = d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3              ]\
			     + d_c[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ];
			FL2y = d_b[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			     + d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4              ];
			FR2y = d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4              ]\
			     + d_c[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ];
			FL3y = d_b[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			     + d_c[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5              ];
			FR3y = d_b[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5              ]\
			     + d_c[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ];

			d_a[ i                              ]\
				= d_a[ i                              ]\
				- Zx * ( FR1x - FL1x ) - Zy * ( FR1y - FL1y );
			d_a[ i + ( Nx + 2 ) * ( Ny + 2)     ]\
				= d_a[ i + ( Nx + 2 ) * ( Ny + 2)     ]\
				- Zx * ( FR2x - FL2x ) - Zy * ( FR2y - FL2y );
			d_a[ i + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
				= d_a[ i + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
				- Zx * ( FR3x - FL3x ) - Zy * ( FR3y - FL3y );	
	}	
	//wait for a while
	__syncthreads();
	//set reflective boundary
	if( i % ( Nx + 2 ) == 0 ){
		d_a[ i                               ] \
		=   d_a[ i + 1                               ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		= - d_a[ i + 1 + ( Nx + 2 ) * ( Ny + 2 )     ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		=   d_a[ i + 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i % ( Nx + 2 ) == Nx + 1 ){
		d_a[ i                               ] \
		=   d_a[ i - 1                               ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		= - d_a[ i - 1 + ( Nx + 2 ) * ( Ny + 2 )     ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		=   d_a[ i - 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i / ( Nx + 2 ) == 0 ){
		d_a[ i                               ] \
		=   d_a[ i + ( Nx + 2 )                               ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		=   d_a[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 )     ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		= - d_a[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i / ( Nx + 2 ) == Ny + 1 ){
		d_a[ i                               ] \
		=   d_a[ i - ( Nx + 2 )                               ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		=   d_a[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 )     ];
		d_a[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		= - d_a[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
}

void Sent_To_Host( float *h_a , float *d_a ){
	cudaError_t Error;
	size_t size;
	printf( "Sending data to Host(CPU)....\n" );
	size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
	Error = cudaMemcpy( h_a , d_a , size , cudaMemcpyDeviceToHost );
	if(DEBUG) printf( "CUDA Error Memcpy d_a->h_a = %s\n" , cudaGetErrorString( Error ) );
}

void Free( float **h_a , float **d_a , float **d_b , float **d_c ){
	cudaError_t Error;
	printf( "Free memory....\n" );
	free( *h_a );
	Error = cudaFree( *d_a );
	if(DEBUG) printf( "CUDA Error cudaFree(d_a) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaFree( *d_b );
	if(DEBUG) printf( "CUDA Error cudaFree(d_b) = %s\n" , cudaGetErrorString( Error ) );
	Error = cudaFree( *d_c );
	if(DEBUG) printf( "CUDA Error cudaFree(d_c) = %s\n" , cudaGetErrorString( Error ) );
}
