#include "main.h"

void Initiate( float *a , float *b , float *c, float *d);
void Save_Result( float *a, float *b, float *c, float *d);

int main(){
	float *h_u, *h_rho, *h_v, *h_T;
	float *d_u, *d_rho, *d_v, *d_T, *d_F, *d_PL, *d_PR, *d_ML, *d_MR;
	int i;
	Allocate_Memory( &h_u, &h_rho, &h_v, &h_T\
		       , &d_u, &d_rho, &d_v, &d_T, &d_F, &d_PL, &d_PR, &d_ML, &d_MR );
	Initiate( h_u, h_rho, h_v, h_T );
	Send_To_Device( h_u, h_rho, h_v, h_T\
		      , d_u, d_rho, d_v, d_T );
	for( i = 0 ; i < NO_STEP ; ++i ){
		GPU_Compute( d_u, d_rho, d_v, d_T, d_F, d_PL, d_PR, d_ML, d_MR );
	}
	Send_To_Host( h_u, h_rho, h_v, h_T\
		    , d_u, d_rho, d_v, d_T );
	Save_Result( h_u, h_rho, h_v, h_T );
	Free( &h_u, &h_rho, &h_v, &h_T\
	    , &d_u, &d_rho, &d_v, &d_T, &d_F, &d_PL, &d_PR, &d_ML, &d_MR );
	return 8;
}

void Initiate( float *a, float *b, float*c, float *d ){
	int i;
	for( i = 0 ; i < N ; ++i ){
		if( i < 0.5 * N ){
			b[ i ] = 10.0;
			c[ i ] = 0.0;
			d[ i ] = 1.0;
		}else{
			b[ i ] = 1.0;
			c[ i ] = 0.0;
			d[ i ] = 1.0;
		}
		a[ i         ] = b[ i ];
		a[ i + N     ] = b[ i ] * c[ i ];
		a[ i + N * 2 ] = b[ i ] * ( CV * d[ i ] + 0.5 * c[ i ] * c[ i ] );
	}
}

void Save_Result( float *a , float *b, float *c, float *d ){
	FILE *output;
	int i;
	output = fopen( "data.txt", "w" );
	if( output == NULL ){
		printf( "File open failed\n" );
	}else{
		for( i = 0; i < N; ++i ){
			fprintf( output, "%e %e %e %e\n", i * DX , b[ i ] , c[ i ], d[ i ] );
		}
	}
	fclose( output );
}
