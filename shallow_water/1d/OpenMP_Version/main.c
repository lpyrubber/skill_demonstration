#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define L       100.0
#define N       20000
#define NT      8
#define DX      (L/N)
#define DT      (0.01*DX)
#define Z       (DT/DX)
#define NO_STEP 80000
#define G       9.81
#define DEBUG   1

void Allocate_Memory( float **a , float **b , float **c , float **d );
void Initial( int i , float *a );
void Compute( float *a , float *b , float *c , float *d );
void Save_Result( float *a );
void Free( float **a , float **b , float **c , float **d );

int main(){
	int np;
	float *u , *u_new , *fm , *fp; 
	double total_t, start_t, end_t;
	int Np , i , j , tid ;

	omp_set_dynamic( 0 );
	omp_set_num_threads( NT );
	start_t=omp_get_wtime();
	tid = omp_get_thread_num();

	Allocate_Memory( &u , &u_new , &fm , &fp );
	#pragma omp parallel shared( u , u_new , fm , fp ) \
			     private( Np , i , j , tid )
	{
		tid = omp_get_thread_num();
		Initial( i , u );	
		#pragma omp barrier
		for ( i = 0 ; i < NO_STEP ; ++i ){
			Compute( u , u_new , fm , fp ); 			
		}
		
	}
	#pragma omp master
	{
		Save_Result( u );
	}
	Free( &u , &u_new , &fm , &fp );
	end_t=omp_get_wtime();
        total_t = end_t-start_t;
        printf("CPU runtime = %lf sec\n", total_t);
	return  0;
}

void Allocate_Memory( float **a , float **b , float **c , float **d ){
	size_t size = 2 * ( N + 2 ) * sizeof( float );
	*a = ( float* )malloc( size );
	*b = ( float* )malloc( size );
	*c = ( float* )malloc( size );
	*d = ( float* )malloc( size );
}

void Initial( int i , float *a ){
	#pragma omp for
	for( i = 0 ; i < N + 2 ; ++i ){
		if( ( i - 1 ) < 0.5 * N ) {
			a[ i         ] = 10.0;
			a[ i + N + 2 ] = 0;
		}else{
			a[ i         ] = 1.0;
			a[ i + N + 2 ] = 0;
		}
	}
}

void Compute( float *a , float *b , float *c , float *d ){
	int i, tid;
	float vel , acc , F1 , F2 , Fr , FL1 , FL2 , FR1 , FR2 ;
	tid = omp_get_thread_num();
	#pragma omp for
	for( i = 0 ; i < N + 2 ; ++i ){
		vel = a[ i + N + 2 ] / a[ i ] ;
		acc = sqrt( G * a[ i ] );
		Fr = vel / acc;
		if ( Fr >   1 ) Fr = 1;
		if ( Fr < - 1 ) Fr = -1 ;
		F1 = a[ i ] * vel ;
		F2 = a[ i ] * vel * vel + 0.5 * G * a[ i ] * a[ i ];
		c[ i         ] =  0.5 * ( F1 * ( Fr + 1 ) + a[ i         ] * acc * ( 1 - Fr * Fr ) ) ;
                c[ i + N + 2 ] =  0.5 * ( F2 * ( Fr + 1 ) + a[ i + N + 2 ] * acc * ( 1 - Fr * Fr ) ) ;
                d[ i         ] = -0.5 * ( F1 * ( Fr - 1 ) + a[ i         ] * acc * ( 1 - Fr * Fr ) ) ;
                d[ i + N + 2 ] = -0.5 * ( F2 * ( Fr - 1 ) + a[ i + N + 2 ] * acc * ( 1 - Fr * Fr ) ) ;
//		if( DEBUG ){
//			printf( "%d %lf %lf %lf %lf %lf\n" , i , Fr , F1 , F2 , vel , acc );
//		}
	}
	
	#pragma omp barrier
	#pragma omp for
	for( i = 1 ; i < N + 1 ; ++i ){
		FL1 = c[ i - 1 ] + d[ i     ];
                FR1 = c[ i     ] + d[ i + 1 ];
                FL2 = c[ i + N + 2 - 1 ] + d[ i + N + 2     ];
                FR2 = c[ i + N + 2     ] + d[ i + N + 2 + 1 ];
                b[ i         ] = a[ i         ] - Z * ( FR1 - FL1 ) ;
                b[ i + N + 2 ] = a[ i + N + 2 ] - Z * ( FR2 - FL2 ) ;
	}
	#pragma omp barrier
	#pragma omp master
	{
		//boundary condition
		b[ 0     ] =   b[         1 ];
		b[ N + 2 ] = - b[ N + 2 + 1 ];
		b[ N         + 1 ] =   b[ N         ];
		b[ N + N + 2 + 1 ] = - b[ N + N + 2 ];
	}
	#pragma omp barrier
	#pragma omp for
	for( i = 0 ; i < N + 2 ; ++i ){
		a[ i         ] = b[ i         ] ;
		a[ i + N + 2 ] = b[ i + N + 2 ] ;
	}
}

void Save_Result( float *a ){
	FILE *pFile;
	int i;
	pFile = fopen( "data_0.dat" , "w" ) ;
	for( i = 1 ; i < N + 1 ; ++i ){
		fprintf( pFile , "%lf %lf %lf\n", i*DX , a[ i ] , a[ i + N + 2 ] );
	}
	fprintf( pFile , "\n" );
	fclose( pFile );
}

void Free( float **a , float **b , float **c , float **d ){
	free(*a);
	free(*b);
	free(*c);
	free(*d);
}
