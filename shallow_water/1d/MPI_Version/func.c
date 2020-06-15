#include "main.h"

void Allocate_Memory( int id , int max , float **a , float **b , float **c ){
	//define variable in function
	int i , averow , extra , rows ;
	if( id == MASTER ){
		*a = ( float* )malloc( N * sizeof( float ) );
		*b = ( float* )malloc( N * sizeof( float ) );
		*c = ( float* )malloc( 1 * sizeof( float ) );
	}else{
		averow = N / ( max - 1 );
		extra  = N % ( max - 1 );
		rows = ( id <= extra ) ?  averow + 3 : averow + 2 ; // +2 for ghost blocks
		*a = ( float* )malloc( rows * sizeof( float ) );
		*b = ( float* )malloc( rows * sizeof( float ) );
		*c = ( float* )malloc( rows * sizeof( float ) );	
	}
}

void Compute( int id , int ele , float *a , float *b ){

}

void Free_Memory( float **a , float **b , float **c ){
	free( *a );
	free( *b );
	free( *c );
}
