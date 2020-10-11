#include "main.h"
#include "mpi.h"

void Initial( float *a );
void Send_To_Slave( int id , int max , int *ele , int *ir , int *il \
		, float *a );
void Communicate( int id , int max , int ele , int ir , int il \
		, float *a );
void Send_To_Master( int id , int max, int ele , float *a );
void Save_Result( float *a );

int main( int argc, char *argv[]){
	//Degin MPI
	MPI_Init( &argc , &argv );
	//set up scalar variable
	int taskid , numtask , left , right , offset , rows;
	int i;
	float *u , *fp , *fm;
	//get variable value from MPI
	MPI_Comm_size( MPI_COMM_WORLD , &numtask );
	MPI_Comm_rank( MPI_COMM_WORLD , &taskid  );

	//allocate Memory
	Allocate_Memory( taskid , numtask , &u , &fp , &fm );
	//only master do the initialization
	if( taskid == MASTER ) Initial( u );
	//distribute  value to each node
	Send_To_Slave( taskid , numtask , &rows , &right , &left \
		     , u );
	for( i = 0 ; i < NO_STEP ; ++i ){
		//exchange boundary between each node for each step
		Communicate( taskid , numtask , rows, right , left \
			   , u );
		Compute( taskid , rows , u , fp , fm );
		Communicate( taskid , numtask , rows, right , left \
			   , u );
	}
	//send all data back to master
	Send_To_Master( taskid , numtask , rows , u );
	if ( taskid == MASTER ) Save_Result(u);
	Free_Memory( &u , &fp , &fm );
	MPI_Finalize();
}

void Initial( float *a ){
	int i , j;
	for( j = 0 ; j < Ny + 2 ; ++j ){
		for( i = 0 ; i < Nx + 2 ; ++i){
			if( ( ( i - 1 ) < 0.5 * Nx ) && ( ( j - 1 ) < 0.5 * Ny ) ){
				a[ i + ( Nx + 2 ) * j                               ] = 10.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] = 0.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] = 0.0;
			}else{
				a[ i + ( Nx + 2 ) * j                               ] = 1.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] = 0.0;
				a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] = 0.0;
			}

/*			a[ i + ( Nx + 2 ) * j                               ] = j;
			a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] = 0.0;
			a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] = 0.0;
*/		}
	}
}

void Send_To_Slave( int id , int max , int *ele , int *ir , int *il\
		, float *a ){
	int offset , rows   , left   , right , dest \
	    ,  tag , source , averow , extra ;
	tag    = BEGIN;
	source = MASTER;
	averow = Ny / (max - 1 );
	extra  = Ny % (max - 1 );
	offset = Nx + 2;
	if( id == MASTER ){
		*ir = NONE;
		*il = NONE;
		//send data to slaves
		for ( dest = 1 ; dest < max ; ++dest ){
			rows = ( dest <= extra ) ? ( averow + 1 ) * ( Nx + 2 ) : averow * ( Nx + 2 );
			//start to send data to slaves
			MPI_Send( a + offset                               , rows , MPI_FLOAT , dest \
				, tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ( Ny + 2 ) + offset     , rows , MPI_FLOAT , dest \
				, tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ( Ny + 2 ) * 2 + offset , rows , MPI_FLOAT , dest \
				, tag , MPI_COMM_WORLD );
			offset += rows;
		}
	}else{
		*il = ( id == 1           ) ? NONE : id - 1;
		*ir = ( id == ( max - 1 ) ) ? NONE : id + 1;
		*ele = ( id <= extra ) ? ( averow + 1 ) : averow ;
		MPI_Recv( a + offset                              , ( *ele ) * ( Nx + 2 ) , MPI_FLOAT , source \
			, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
		MPI_Recv( a + ( Nx + 2 ) * ( Ny + 2 ) + offset     , ( *ele ) * ( Nx + 2 ) , MPI_FLOAT , source \
			, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
		MPI_Recv( a + ( Nx + 2 ) * ( Ny + 2 ) * 2 + offset , ( *ele ) * ( Nx + 2 ) , MPI_FLOAT , source \
			, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
	}
}

void Send_To_Master( int id , int max , int ele , float *a ){
	int tag , dest , source , extra , averow , offset , rows;
	tag    = FINAL;
	dest   = MASTER;
	averow = Ny / ( max - 1 );
	extra  = Ny % ( max - 1 );
	offset = Nx + 2;
	if( id == MASTER ){
		for( source = 1 ; source < max ; ++source ){
			rows = ( source <= extra ) ? ( averow + 1 ) * ( Nx + 2 ) : averow * ( Nx + 2 );
			MPI_Recv( a + offset                               , rows , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( Ny + 2 ) + offset     , rows , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( Ny + 2 ) * 2 + offset , rows , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			offset += rows;
		}
	}else{
		MPI_Send( a + offset                                , ele * (Nx + 2 ) , MPI_FLOAT , dest , tag \
			, MPI_COMM_WORLD );
		MPI_Send( a + ( Nx + 2 ) * ( ele + 2 )     + offset , ele * (Nx + 2 ) , MPI_FLOAT , dest , tag \
			, MPI_COMM_WORLD );
		MPI_Send( a + ( Nx + 2 ) * ( ele + 2 ) * 2 + offset , ele * (Nx + 2 ) , MPI_FLOAT , dest , tag \
			, MPI_COMM_WORLD );
	}
}

void Communicate( int id   , int max  , int ele , int ir , int il 
		, float *a ){
	//set variables inside function
	int tag , source , dest , i , j;
	if( id != MASTER ){
		//excchange between odd & even
		//  ______________________
		//  |  |  |  |  |  |  |  |
		//  |__|__|__|__|__|__|__|
		//   1  2  3  4  5  6  7
		//    <->   <->   <->   <x>
		//        even id    left core exist
		if ( ( id % 2 == 0 ) && ( il != NONE ) ){
			tag = LEFT ; //send->recv in left direction
			dest = il ;
			//send boundary to left core
			MPI_Send( a + Nx + 2                               , Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
			MPI_Send( a + Nx + 2 + ( Nx + 2 ) * ( ele + 2 )    , Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
			MPI_Send( a + Nx + 2 + ( Nx + 2 ) * ( ele + 2 ) * 2, Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
			tag = RIGHT ;
			source = il ; //send->recv in right direction
			//recv boundary from left core
			MPI_Recv( a                                , Nx + 2 , MPI_FLOAT , source \
				, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 , MPI_FLOAT , source \
				, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 , MPI_FLOAT , source \
				, tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
		}
		//        odd id     right core exist      
		if ( ( id % 2 == 1 ) && ( ir != NONE ) ){
			tag = LEFT ; //send->recv in left direction
			source = ir ; 
			//receive boundary from right core
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 )                                , Nx + 2 \
				, MPI_FLOAT , source   , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 \
				, MPI_FLOAT , source   , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 \
				, MPI_FLOAT , source   , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			tag = RIGHT ; //exchange in right direction
			dest = ir ;
			//send boundary to right core
			MPI_Send( a + ( Nx + 2 ) * ele                                 , Nx + 2 , MPI_FLOAT \
				, dest , tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ele  + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 , MPI_FLOAT \
				, dest , tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ele  + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 , MPI_FLOAT \
				, dest , tag , MPI_COMM_WORLD );
		}

		//excchange between even & odd
		//  ______________________
		//  |  |  |  |  |  |  |  |
		//  |__|__|__|__|__|__|__|
		//   1  2  3  4  5  6  7
		// <x>   <->   <->   <->   
		//       even id      right core exist
		if ( ( id % 2 == 0 ) && ( ir != 0 ) ){
			tag = RIGHT ; //send->recv in right direction
			dest = ir ;
			//send boundary to left core
			MPI_Send( a + ( Nx + 2 ) * ele                               , Nx + 2 , MPI_FLOAT \
				, dest   , tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ele + ( Nx + 2 ) * ( ele + 2 )    , Nx + 2 , MPI_FLOAT \
				, dest   , tag , MPI_COMM_WORLD );
			MPI_Send( a + ( Nx + 2 ) * ele + ( Nx + 2 ) * ( ele + 2 ) * 2, Nx + 2 , MPI_FLOAT \
				, dest   , tag , MPI_COMM_WORLD );
			tag = LEFT ; //send->recv in left direction
			source = ir ;
			//recv boundary from left core
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 )                                , Nx + 2 , MPI_FLOAT \
				, source , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 , MPI_FLOAT \
				, source , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 , MPI_FLOAT \
				, source , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE );
		}
		//        odd  id       lefr core exist      
		if ( ( id % 2 == 1 ) && ( il != 0 ) ){
			tag = RIGHT ; //send->recv in right direction
			source = il ;
			//receive boundary from right core
			MPI_Recv( a                                , Nx + 2 , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			MPI_Recv( a + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 , MPI_FLOAT , source , tag \
				, MPI_COMM_WORLD , MPI_STATUS_IGNORE );
			tag = LEFT ; //send->recv in left direction
			dest = il ;
			//send boundary to right core
			MPI_Send( a + Nx + 2                                , Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
			MPI_Send( a + Nx + 2 + ( Nx + 2 ) * ( ele + 2 )     , Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
			MPI_Send( a + Nx + 2 + ( Nx + 2 ) * ( ele + 2 ) * 2 , Nx + 2 , MPI_FLOAT , dest   , tag \
				, MPI_COMM_WORLD );
		}
		//set reflective boundary
		if( id == 1 ){
			for( i = 0 ; i < Nx + 2 ; ++i ){
				a[ i                                ] = \
				  a[ i + Nx + 2                               ];
				a[ i + ( Nx + 2 ) * ( ele + 2 )     ] = \
				  a[ i + Nx + 2 + ( Nx + 2 ) * ( ele + 2 )    ];
				a[ i + ( Nx + 2 ) * ( ele + 2 ) * 2 ] = \
				- a[ i + Nx + 2 + ( Nx + 2 ) * ( ele + 2 ) * 2];
			}
		}else if( id == max - 1 ){
			for( i = 0 ; i < Nx + 2 ; ++i ){
				a[ i + ( Nx + 2 ) * ( ele + 1 )                                ] = \
				  a[ i + ( Nx + 2 ) * ele                                ];
				a[ i + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 )     ] = \
				  a[ i + ( Nx + 2 ) * ele + ( Nx + 2 ) * ( ele + 2 )     ];
				a[ i + ( Nx + 2 ) * ( ele + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 2 ] = \
				- a[ i + ( Nx + 2 ) * ele + ( Nx + 2 ) * ( ele + 2 ) * 2 ];
			}
		}
		for ( j = 0 ; j < ele + 2 ; ++j ){
			a[ ( Nx + 2 ) * j                                ] = \
			  a[ 1 + ( Nx + 2 ) * j                                ];
			a[ ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] = \
			- a[ 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ];
			a[ ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] = \
			  a[ 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ];

			a[ Nx + 1 + ( Nx + 2 ) * j                                ] = \
			  a[ Nx + ( Nx + 2 ) * j                                ];
			a[ Nx + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] = \
			- a[ Nx + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ];
			a[ Nx + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] = \
			  a[ Nx + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ];			
		}

	}
}

void Save_Result( float *a ){
	int i , j;
	FILE *pFile;
        
        pFile = fopen( "data.txt", "w" );
	for( j = 1 ; j < Ny + 1 ; ++j){
 	       	for( i = 1 ; i < Nx + 1 ; ++i){
	       		fprintf(pFile, "%e %e %e %e %e\n"\
				, i * DX\
				, j * DY\
				, a[ i + ( Nx + 2 ) * j                               ]\
				, a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
				, a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]);
                }
	}
	fclose(pFile);
}
