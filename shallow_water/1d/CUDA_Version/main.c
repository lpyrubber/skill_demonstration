#include "main.h"

void Initiate( float *a , float *b );
void Save_Result( float *a , float *b );

int main(){
	//set variable
	float *h_u , *h_s , *d_u , *d_s , *d_f ;
	char  buffer[ 20 ];
	int  i;
        double total_t;
        clock_t start_t, end_t;
        start_t=clock();
	Allocate_Memory( &h_u , &h_s , &d_u , &d_s , &d_f );
	Initiate( h_u , h_s );
	Sent_To_Device( h_u , h_s , d_u , d_s );
	for( i = 0 ; i < NO_STEP ; ++i ){
		GPU_Compute( d_u , d_s , d_f );
	}
	Sent_To_Host( h_u , h_s , d_u , d_s );
	Save_Result( h_u , h_s );
	Free( &h_u , &h_s , &d_u , &d_s , &d_f );
	end_t=clock();
        total_t=(double)(end_t-start_t)/CLOCKS_PER_SEC;
        printf("Total time taken by CPU: %lf sec\n", total_t);
	return 0;
}

void Initiate( float *a , float *b ){
	int i;
	for( i = 0 ; i < N + 2 ; ++i ){
 		if( ( i - 1 ) < 0.5 * N ){
			a[ i ] = 10.0;
			b[ i ] =  0.0;
		}else{
			a[ i ] =  1.0;
			b[ i ] =  0.0;
		}
	}
}

void Save_Result( float *a , float *b ){
	FILE *output;
	int i;
	output = fopen( "data_0.dat", "w" );
	if( output == NULL ){
		printf( "File open failed\n" );
	}else{
		for( i = 1 ; i < N + 1 ; ++i ){
			fprintf( output , "%e %e %e\n" , i*DX , a[i] , b[i] );
		}
	}
	fclose( output );
}


