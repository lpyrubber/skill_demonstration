#include "main.h"
#include "mpi.h"

void Initial( float *a , float *b , float *c , float *d );
void Sent_To_Slave( int   id , int  max , int *ele , int *ir , int *il \
		  , float *a , float *b , float *c , float *d);
void Communicate( int id   , int max , int ele , int ir , int il \
		, float *a , float *b , float *c , float *d);
void Sent_To_Master( int id   , int max , int ele, float *a , float *b , float *c , float *d );
void Save_Result( float *a , float *b , float *c);

int main( int argc , char *argv[] ){
	//Begin MPI
	MPI_Init( &argc , &argv );	
	//note : each core while only read its owe variable
	//set up scalar variable
	int taskid , numtask , left , right , offset , rows \
	  , i      , j       , k;
	//set up array variable by pointer
	float *u , *rho , *v , *T , *f , *PL , *PR , *ML , *MR;
	//get variable value from MPI
	MPI_Comm_size( MPI_COMM_WORLD , &numtask );
	MPI_Comm_rank( MPI_COMM_WORLD , &taskid );
	//print debug message
	if( DEBUG ) printf( "num of task = %d , task id = %d\n" , numtask , taskid );
	//allocate memory
	Allocate_Memory( taskid , numtask , &u , &rho , &v , &T , &f , &PL , &PR , &ML , &MR );
	//only master core do the initialization
	if( taskid == MASTER ) Initial(u, rho , v , T );
	//distribute value to each core
	Sent_To_Slave( taskid , numtask , &rows, &right , &left \
		     , u , rho , v , T );
//	if(DEBUG) printf( "id = %d, right = %d , left = %d\n" , taskid, right, left );
	//start to compute
	for( i = 0 ;  i < NO_STEP ; ++i ){
		//exchange boundary value between each core for each step
		Communicate( taskid , numtask , rows , right , left \
			   , u , rho , v , T );
		Compute( taskid , rows , u , rho , v , T , f , PL , PR , ML , MR );
	}
	//sent all the value back to master
	Sent_To_Master( taskid , numtask , rows , u , rho , v , T );
	//master save & print the result
	if( taskid == MASTER ) Save_Result(rho, v , T );
	//clear the dynamic memory
	Free_Memory( &u , &rho , &v , &T , &f , &PL , &PR , &ML , &MR );
	//close MPI
	MPI_Finalize();
}

//initiate process
void Initial( float *a, float *b , float *c , float *d ){
	int i;
	for( i = 0 ; i < N+1 ; ++i){ //plus two for the ghost blocks
		if( ( ( i - 1 ) < 0.5 * N ) ){
			b[ i ] = 10.0; //rho
			c[ i ] = 0;    //rho * u
			d[ i ] = 1.0;
		}else{
			b[ i ] = 1.0;  //rho
			c[ i ] = 0;    //rho * u
			d[ i ] = 1.0;
		}
		a[ i ] = b[ i ];
		a[ i + N + 1 ] = b[ i ] * c[ i ];
		a[ i + ( N + 1 ) * 2 ] = b[ i ] * ( CV * d[ i ] + 0.5 * c[ i ] * c[ i ] );
	}
}

void Sent_To_Slave( int   id , int   max , int *ele , int *ir , int *il \
		  , float *a , float *b , float *c , float *d){
	int offset , rows   , left   , right , dest \
	  , tag    , source , averow , extra ;
	//using tag to identify which pair of send & recv for the same dest
	// & source
	tag = BEGIN ;
	//compute basic value
	averow = N / max;
	extra  = N % max;
	offset = 0;
	//there is nothing left for the leftest core
	*il  = ( id == MASTER      ) ? NONE : id - 1 ;         
	//there is nothing right for the rightest core
	*ir  = ( id == ( max - 1 ) ) ? NONE : id + 1 ;    
	*ele = ( id < extra ) ? averow + 1 : averow ;   
	if( id == MASTER ){ //MASTER is fixed to be core 0
		offset = *ele + 1;
		//send data to slaves
		for ( dest = 1 ; dest < max ; ++dest ){
			rows = ( dest < extra ) ? averow + 1 : averow ;
			//start to send data to slaves
			MPI_Send( a + offset                 , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( a + offset + N + 1         , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( a + offset + ( N + 1 ) * 2 , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( b + offset , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( c + offset , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			MPI_Send( d + offset , rows , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
			offset += rows ;
			//show debug message
//			if( DEBUG ) printf( "send data to %d\n" , dest );	
		}	
		*ele = ( id < extra ) ? averow + 1 : averow ;   
	}else{
		source = MASTER;
		//start to receive data from master
                MPI_Recv( a + 1                  , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + 1 + *ele + 2         , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + 1 + ( *ele + 2 ) * 2 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
                MPI_Recv( b + 1 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
                MPI_Recv( c + 1 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
                MPI_Recv( d + 1 , *ele , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );

		//show debug message
//		if( DEBUG ) printf( "%d recv data\n" , id );
	}	
}

void Sent_To_Master( int id , int max , int ele , float *a, float *b , float *c , float *d ){
	//define the variable inside this function
	int tag , dest , source , extra , averow , offset , rows ; 
	tag = FINAL ;
	dest = MASTER;
	if( id == MASTER ){
		//master count rows & offset for each core
		averow = N / max ;
        	extra  = N % max ;
		rows = ( source < extra ) ? averow + 1 : averow ;
       		offset = rows + 1;
		for( source = 1 ; source < max ; ++source){
			rows = ( source < extra ) ? averow + 1 : averow ;
			MPI_Recv( a + offset                 , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( a + N + 1 + offset         , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( a + ( N + 1 ) * 2 + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( b + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( c + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			MPI_Recv( d + offset , rows , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
				, MPI_STATUS_IGNORE);
			offset += rows;
			if(DEBUG) printf( "receive from %d, %d ,%d\n" , source , rows, offset); 
		}
	}else{
		//slaves already have their own information, therefore, they just send 
		//array to master
		MPI_Send( a + 1                   , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( a + 1 + ele + 2         , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( a + 1 + ( ele + 2 ) * 2 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( b + 1 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( c + 1 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);
		MPI_Send( d + 1 , ele , MPI_FLOAT , dest , tag , MPI_COMM_WORLD);

//		if(DEBUG) printf( "%d send to master\n" , id );
	}
}

void Communicate( int id   , int max  , int ele , int ir , int il 
		, float *a , float *b , float *c , float *d ){
	//set variables inside function
	int tag , source , dest , len ;
	len = ( id == MASTER ) ? N + 1 : ele + 2;
	//excchange between odd & even
	//  _________________________
	//  |  |  |  |  |  |  |  |  |
	//  |__|__|__|__|__|__|__|__|
	//   0  1  2  3  4  5  6  7
	//    <->   <->   <->   <x>
	//        odd id     left core exist
	if ( ( id % 2 == 1 ) && ( il != NONE ) ){
		tag = LEFT ; //send->recv in left direction
		dest = il ;
		//send boundary to left core
		MPI_Send( a + 1           , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + 1 + len     , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + 1 + len * 2 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( b + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( c + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( d + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		tag = RIGHT ;
		source = il ; //send->recv in right direction
		//recv boundary from left core
		MPI_Recv( a           , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + len     , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + len * 2 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( b , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( c , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( d , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );

	}
	//        even id    right core exist      
	if ( ( id % 2 == 0 ) && ( ir != NONE ) ){
		tag = LEFT ; //send->recv in left direction
		source = ir ; 
		//receive boundary from right core
		MPI_Recv( a + ele + 1          , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		MPI_Recv( a + ele + 1 + len    , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		MPI_Recv( a + ele + 1 + len * 2, 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		MPI_Recv( b + ele + 1 , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		MPI_Recv( c + ele + 1 , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		MPI_Recv( d + ele + 1 , 1 , MPI_FLOAT , source   , tag , MPI_COMM_WORLD \
   				, MPI_STATUS_IGNORE );
		tag = RIGHT ; //exchange in right direction
		dest = ir ;
		//send boundary to right core
		MPI_Send( a + ele           , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
		MPI_Send( a + ele + len     , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
		MPI_Send( a + ele + len * 2 , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
		MPI_Send( b + ele , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
		MPI_Send( c + ele , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
		MPI_Send( d + ele , 1 , MPI_FLOAT , dest , tag , MPI_COMM_WORLD );
	}

	//excchange between even & odd
	//  _________________________
	//  |  |  |  |  |  |  |  |  |
	//  |__|__|__|__|__|__|__|__|
	//   0  1  2  3  4  5  6  7
	// <x>   <->   <->   <->   
	//       odd id      right core exist
	if ( ( id % 2 == 1 ) && ( ir != NONE ) ){
		tag = RIGHT ; //send->recv in right direction
		dest = ir ;
		//send boundary to left core
		MPI_Send( a + ele           , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + ele + len     , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + ele + len * 2 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( b + ele , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( c + ele , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( d + ele , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		tag = LEFT ; //send->recv in left direction
		source = ir ;
		//recv boundary from left core
		MPI_Recv( a + ele + 1 , 1          , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + ele + 1 + len    , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + ele + 1 + len * 2, 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( b + ele + 1 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( c + ele + 1 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( d + ele + 1 , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
	}
	//        even  id       lefr core exist      
	if ( ( id % 2 == 0 ) && ( il != NONE ) ){
		tag = RIGHT ; //send->recv in right direction
		source = il ;
		//receive boundary from right core
		MPI_Recv( a          , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + len    , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( a + len * 2, 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( b , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( c , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		MPI_Recv( d , 1 , MPI_FLOAT , source , tag , MPI_COMM_WORLD \
			, MPI_STATUS_IGNORE );
		tag = LEFT ; //send->recv in left direction
		dest = il ;
		//send boundary to right core
		MPI_Send( a + 1           , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + 1 + len     , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( a + 1 + len * 2 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( b + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( c + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
		MPI_Send( d + 1 , 1 , MPI_FLOAT , dest   , tag , MPI_COMM_WORLD );
	}
	//set reflective boundary
	if( id == 0 ){
		a[ 0       ] =  a[ 1           ];
		a[ len     ] = -a[ 1 + len     ];
		a[ len * 2 ] =  a[ 1 + len * 2 ];
		b[ 0 ] =   b[ 1 ];
		c[ 0 ] = - c[ 1 ];
		d[ 0 ] =   d[ 1 ];
	}else if( id == max - 1 ){
		a[ ele + 1          ] =  a[ ele           ];
		a[ ele + 1 + len    ] = -a[ ele + len     ];
		a[ ele + 1 + len * 2] =  a[ ele + len * 2 ];
		b[ ele + 1 ] =   b[ ele ];
		c[ ele + 1 ] = - c[ ele ];
		d[ ele + 1 ] =   d[ ele ];
	}
}

void Save_Result( float *a , float *b , float *c){
	int i ;
	static int time = 0 ;
	FILE *output ;
	char buffer[ 15 ];
	
	//create file name by static integral
	sprintf( buffer , "data.txt" );
	//open output file
	output = fopen( buffer , "w");
	for( i = 1 ; i < N+1 ; ++i ){
	       fprintf( output , "%lf %lf %lf %lf\n" , DX * i , a[ i ] , b[ i ] , c[ i ]);
	}
	fclose ( output );
	++time;
}
