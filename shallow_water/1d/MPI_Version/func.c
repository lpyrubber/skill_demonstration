#include "main.h"

void Allocate_Memory( int id , int max , float **a , float **b , float **c ){
	//define variable in function
	int i , averow , extra , rows ;
	averow = N / max;
	extra  = N % max;
	if( id == MASTER ){
		*a = ( float* )malloc( (N+1) * sizeof( float ) );
		*b = ( float* )malloc( (N+1) * sizeof( float ) );
		*c = ( float* )malloc( 4 * ( averow + 3 ) * sizeof( float ) );
	}else{
		rows = ( id < extra ) ?  averow + 3 : averow + 2 ; // +2 for ghost blocks
		*a = ( float* )malloc( rows * sizeof( float ) );
		*b = ( float* )malloc( rows * sizeof( float ) );
		*c = ( float* )malloc( 4 * rows * sizeof( float ) );	
	}
}

void Compute( int id , int ele , float *a , float *b , float *c ){
	//create variable in this function
	int i;
	float vel , acc , F1 , F2 , Fr , FL1 , FL2 , FR1 , FR2 ;
	//compute flux
	for( i = 0 ; i < ele + 2 ; ++i ){
		vel = b[ i ] / a[ i ] ;
		acc = sqrt( G * a[ i ] ) ;
		Fr  = vel / acc ;
		if( Fr >  1 ) Fr =  1 ;
		if( Fr < -1 ) Fr = -1 ;
		F1 = a[ i ] * vel ;
		F2 = a[ i ] * vel * vel + 0.5 * G * a[ i ] \
		   * a[ i ] ;
		c[ i                   ] =  0.5 * ( F1 * ( Fr + 1 ) + a[ i ] \
					 * acc * ( 1 - Fr * Fr ) ); //fp[i]
		c[ i + ( ele + 2 )     ] =  0.5 * ( F2 * ( Fr + 1 ) + b[ i ] \
					 * acc * ( 1 - Fr * Fr ) ); //fp[i+ele+2]
		c[ i + ( ele + 2 ) * 2 ] = -0.5 * ( F1 * ( Fr - 1 ) + a[ i ] \
					 * acc * ( 1 - Fr * Fr ) ); //fm[i]
		c[ i + ( ele + 2 ) * 3 ] = -0.5 * ( F2 * ( Fr - 1 ) + b[ i ] \
					 * acc * ( 1 - Fr * Fr ) ); //fm[i+ele+2]
	}
	//compute u
	for( i = 1 ; i < ele + 1 ; ++i ){
		FL1 = c[ i - 1                   ] + c[ i     + ( ele + 2 ) * 2 ];
		FR1 = c[ i                       ] + c[ i + 1 + ( ele + 2 ) * 2 ];
		FL2 = c[ i - 1 + ( ele + 2 ) * 1 ] + c[ i     + ( ele + 2 ) * 3 ];
		FR2 = c[ i     + ( ele + 2 ) * 1 ] + c[ i + 1 + ( ele + 2 ) * 3 ];
		a[ i ] = a[ i ] - Z * ( FR1 - FL1 );
		b[ i ] = b[ i ] - Z * ( FR2 - FL2 );
	}
}

void Free_Memory( float **a , float **b , float **c ){
	free( *a );
	free( *b );
	free( *c );
}
