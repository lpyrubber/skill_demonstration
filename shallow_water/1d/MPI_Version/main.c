#include "main.h"
#include "mpi.h"

void Initial( float *a , float *b );
void Sent_To_Slave( int   id , int  max , int *ele , int *ir , int *il \
		  , float *a , float *b );
void Communicate( int id   , int max , int ele , int ir , int il \
		, float *a , float *b );
void Sent_To_Master( int id   , int max , int ele, float *a , float *b );
void Save_Result( float *a , float *b );

int main( int argc , char *argv[] ){
	//Begin MPI
	MPI_Init( &argc , &argv );	
	//note : each core while only read its owe variable
	//set up scalar variable
	int taskid , numtask , left , right , offset , rows \
	  , i      , j       , k;
	//set up array variable by pointer
	float *u , *s , *f;
	//get variable value from MPI
	MPI_Comm_size( MPI_COMM_WORLD , &numtask );
	MPI_Comm_rank( MPI_COMM_WORLD , &taskid );
	//print debug message
//	if( DEBUG ) printf( "num of task = %d , task id = %d\n" , numtask , taskid );
	//allocate memory
	Allocate_Memory( taskid , numtask , &u , &s , &f );
	//only master core do the initialization
	if( taskid == MASTER ) Initial(u, s);
	//distribute value to each core
	Sent_To_Slave( taskid , numtask , &rows, &right , &left \
		     , u      , s );
	if(DEBUG) printf( "id = %d, right = %d , left = %d\n" , taskid, right, left );
	//start to compute
	for( i = 0 ;  i < NO_STEP ; ++i ){
		//exchange boundary value between each core for each step
		Communicate( taskid , numtask , rows , right , left \
			   , u , s );
		Compute( taskid , rows , u , s , f );
	}
	//sent all the value back to master
	Sent_To_Master( taskid , numtask , rows , u , s );
	//master save & print the result
	if( taskid == MASTER ) Save_Result(u, s);
	//clear the dynamic memory
	Free_Memory( &u , &s, &f );
	//close MPI
	MPI_Finalize();
}

//initiate process
void Initial( float *a, float *b ){
	int i;
	for( i = 0 ; i < N ; ++i){ //plus two for the ghost blocks
		if( ( ( i - 1 ) < 0.5 * N ) ){
			a[ i ] = 10.0; //rho
			b[ i ] = 0;    //rho * u
		}else{
			a[ i ] = 1.0;  //rho
			b[ i ] = 0;    //rho * u
		}
	}
//	if(DEBUG){
//		for( i = 0 ; i < N ; ++i ){
//			printf( "%3.3e %3.3e\n" , a[i] , b[i] );
//		}
//	}

}

void Sent_To_Slave( int   id , int   max , int *ele , int *ir , int *il \
		  , float *a , float *b ){
	int offset , rows   , left   , right , dest \
	  , tag    , source , averow , extra ;
	//using tag to identify which pair of send & recv for the same dest
	// & source
	tag = BEGIN ;
	//compute basic value
	averow = N / ( max - 1 );
	extra  = N % ( max - 1 );
	offset = 0;
	if( id == MASTER ){ //MASTER is fixed to be core 0
		*ir = NONE ;
		*il = NONE ;
		offset = 0 ;
		//send data to slaves
		for ( dest = 1 ; dest < max ; ++dest ){
			rows = ( dest <= extra ) ? averow + 1 : averow ;
			left = ( dest == 1 ) ? 0 : dest - 1 ;
			right = ( dest == ( max - 1 ) ) ? 0 : dest + 1 ;
			//start to send data to slaves
			MPI_Send( a + offset , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( b + offset , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			offset += rows ;
			//show debug message
//			if( DEBUG ) printf( "send data to %d\n" , dest );	
		}
	}else{
		source = MASTER;
		*ele = ( id <= extra ) ? averow + 1 : averow ; 
		//there is nothing left for the leftest core
		*il  = ( id == 1           ) ? NONE : id - 1 ;         
		//there is nothing right for the rightest core
		*ir  = ( id == ( max - 1 ) ) ? NONE : id + 1 ;    
		//start to receive data from master
                MPI_Recv( a + 1 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
                MPI_Recv( b + 1 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		//show debug message
//		if( DEBUG ) printf( "%d recv data\n" , id );
	}	
}

void Sent_To_Master( int id , int max , int ele , float *a, float *b ){
	//define the variable inside this function
	int tag , dest , source , extra , averow , offset , rows ; 
	tag = FINAL ;
	dest = MASTER;
	if( id == MASTER ){
		//master count rows & offset for each core
		averow = N / ( max - 1 );
        	extra  = N % ( max - 1 );
       		offset = 0;
		for( source = 1 ; source < max ; ++source){
			rows = ( source <= extra ) ? averow + 1 : averow ;
			MPI_Recv( a + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( b + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			offset += rows;
			if(DEBUG) printf( "receive from %d, %d ,%d\n" , source , rows, offset); 
		}
	}else{
		//slaves already have their own information, therefore, they just send 
		//array to master
		MPI_Send( a + 1 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( b + 1 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
//		if(DEBUG) printf( "%d send to master\n" , id );
	}
}

void Communicate( int id   , int max  , int ele , int ir , int il 
		, float *a , float *b ){
	//set variables inside function
	int tag , source , dest;
	if( id != MASTER ){
		//excchange between odd & even
		//  ______________________
		//  |  |  |  |  |  |  |  |
		//  |__|__|__|__|__|__|__|
		//   1  2  3  4  5  6  7
		//    <->   <->   <->   <x>
		//        even id    left core exist
		if ( ( id % 2 == 0 ) && ( il != 0 ) ){
			tag = LEFT ; //send->recv in left direction
			dest = il ;
			//send boundary to left core
			MPI_Send( a + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
			MPI_Send( b + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
			tag = RIGHT ;
			source = il ; //send->recv in right direction
			//recv boundary from left core
			MPI_Recv( a     , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
			MPI_Recv( b     , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
		}
		//        odd id     right core exist      
		if ( ( id % 2 == 1 ) && ( ir != 0 ) ){
			tag = LEFT ; //send->recv in left direction
			source = ir ; 
			//receive boundary from right core
			MPI_Recv( a + ele + 1 , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
    				, MPI_STATUS_IGNORE );
			MPI_Recv( b + ele + 1 , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
    				, MPI_STATUS_IGNORE );
			tag = RIGHT ; //exchange in right direction
			dest = ir ;
			//send boundary to right core
			MPI_Send( a + ele     , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( b + ele     , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
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
			MPI_Send( a + ele     , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
			MPI_Send( b + ele     , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
			tag = LEFT ; //send->recv in left direction
			source = ir ;
			//recv boundary from left core
			MPI_Recv( a + ele + 1 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
			MPI_Recv( b + ele + 1 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
		}
		//        odd  id       lefr core exist      
		if ( ( id % 2 == 1 ) && ( il != 0 ) ){
			tag = RIGHT ; //send->recv in right direction
			source = il ;
			//receive boundary from right core
			MPI_Recv( a     , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
			MPI_Recv( b     , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE );
			tag = LEFT ; //send->recv in left direction
			dest = il ;
			//send boundary to right core
			MPI_Send( a + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
			MPI_Send( b + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		}
		//set reflective boundary
		if( id == 1 ){
			a[ 0 ] =   a[ 1 ];
			b[ 0 ] = - b[ 1 ];
		}else if( id == max - 1 ){
			a[ ele + 1 ] =   a[ ele ];
			b[ ele + 1 ] = - b[ ele ];
		}
	}
}

void Save_Result( float *a , float *b ){
	int i ;
	static int time = 0 ;
	FILE *output ;
	char buffer[ 15 ];
	
	//create file name by static integral
	sprintf( buffer , "data_%d.dat" , time );
	//open output file
	output = fopen( buffer , "w");
	for( i = 0 ; i < N ; ++i ){
	       fprintf( output , "%e %e %e\n" , DX * i , a[ i ] , b[ i ]);
	}
	fclose ( output );
	++time;
}

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
		*c = ( float* )malloc( 4 * rows * sizeof( float ) );	
	}
}

void Compute( int id , int ele , float *a , float *b , float *c ){
	//create variable in this function
	int i;
	float vel , acc , F1 , F2 , Fr , FL1 , FL2 , FR1 , FR2 ;
	if( id != MASTER){
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
}

void Free_Memory( float **a , float **b , float **c ){
	free( *a );
	free( *b );
	free( *c );
}
