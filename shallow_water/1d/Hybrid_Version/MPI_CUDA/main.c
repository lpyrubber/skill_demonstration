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
	float *h_u, *h_s ,*d_u , *d_s , *d_f;
	//get variable value from MPI
	MPI_Comm_size( MPI_COMM_WORLD , &numtask );
	MPI_Comm_rank( MPI_COMM_WORLD , &taskid );
	//print debug message
	if( DEBUG ) printf( "num of task = %d , task id = %d\n" , numtask , taskid );
	//allocate memory
	Allocate_Memory( taskid , numtask , &h_u , &h_s , &d_u , &d_s , &d_f );
	//only master core do the initialization
	if( taskid == MASTER ) Initial( h_u, h_s );
	//distribute value to each core
	Sent_To_Slave( taskid , numtask , &rows, &right , &left \
		     , h_u      , h_s );
//	if(DEBUG) printf( "id = %d, right = %d , left = %d\n" , taskid, right, left );
	//start to compute
	for( i = 0 ;  i < NO_STEP ; ++i ){
		//exchange boundary value between each core for each step
		Communicate( taskid , numtask , rows , right , left \
			   , h_u , h_s );
		//send to gpu
		Sent_To_Device( rows , h_u , h_s , d_u , d_s );
		Compute( rows , d_u , d_s , d_f );
		Sent_To_Host( rows , h_u , h_s , d_u, d_s );
		//send it back
	}
	//sent all the value back to master
	Sent_To_Master( taskid , numtask , rows , h_u , h_s );
	//master save & print the result
	if( taskid == MASTER ) Save_Result(h_u, h_s);
	//clear the dynamic memory
	Free_Memory( taskid , &h_u , &h_s , &d_u , &d_s, &d_f );
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
	       fprintf( output , "%lf %lf %lf\n" , DX * i , a[ i ] , b[ i ] / a[ i ]);
	}
	fclose ( output );
	++time;
}
