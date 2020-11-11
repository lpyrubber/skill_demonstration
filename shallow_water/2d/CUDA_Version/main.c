#include "main.h"

void Initiate( float *a);
void Save_Result( float *a);

int main(){
	//set variable
	float *h_u, *d_u , *d_fm , *d_fp ;
	char  buffer[ 20 ];
	int  i ;
	Allocate_Memory( &h_u , &d_u , &d_fm , &d_fp );
	Initiate( h_u );
	Sent_To_Device( h_u , d_u );
	for( i = 0 ; i < NO_STEP ; ++i ){
		GPU_Compute( d_u , d_fm , d_fp );
	}
	Sent_To_Host( h_u , d_u );
	Save_Result( h_u );
	Free( &h_u , &d_u , &d_fm , &d_fp );
	return 0;
}

void Initiate( float *a ){
	int i, j;
	for( j = 0 ; j < Ny + 2 ; ++j ){
		for( i = 0 ; i < Nx + 2 ; ++i ){
 			if( ( j - 1 ) < 0.5 * Ny ){
				a[ i + ( Nx + 2 ) * j                              ] = 10.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)     ] = 0.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2) * 2 ] = 0.0;
			}else{
				a[ i + ( Nx + 2 ) * j                              ] = 1.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)     ] = 0.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2) * 2 ] = 0.0;
			}
		}
	}
}

void Save_Result( float *a ){
	FILE *output;
	int i, j;
	output = fopen( "data_0.dat", "w" );
	if ( output == NULL){
                printf("File open failed\n");
        }else{
		for( j = 1 ; j < Ny + 1 ; ++j){
 	        	for( i = 1 ; i < Nx + 1 ; ++i){
	        		fprintf( output , "%e %e %e %e %e\n"\
					, i * DX\
					, j * DY\
					, a[ i + ( Nx + 2 ) * j                               ]\
					, a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
					, a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]);
        	        }
		}
        }
	fclose( output );
}


