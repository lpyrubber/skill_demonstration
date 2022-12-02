#include "main.h"

int main(){
	float *x, *y, *Fx, *Fy, *K;
	float *Kr, *Fr, *U, *u;
	float *S, scale = 1e4;
	int np, ne, nr, i, j;
	int *tri;
	char *Fix_x, *Fix_y;
	Count_Number( &np, &ne );
	Allocate_Memory( np, ne, &Fix_x, &Fix_y, &tri, &x, &y, &Fx, &Fy, &K, &Kr, &Fr, &u, &U, &S);
	Initial( np, ne, Fix_x, Fix_y, tri, x, y, Fx, Fy);
	Calculate_K(np, ne, tri, x, y, K);
	Partition_K_F( np, &nr , Fix_x, Fix_y , Fx, Fy, Fr , K , Kr);
	Conjugate_Gradient(Kr, Fr, u, nr);
	Partition_U(np, nr, Fix_x , Fix_y , u, U);
	for( i = 0 ; i < np ; ++i ){
		x[i]+=scale*U[2*i];
		y[i]+=scale*U[2*i+1];
	}
	Save_Result( np, ne , tri,  x , y);
	Free_Memory( &Fix_x, &Fix_y, &tri, &x, &y, &Fx, &Fy, &K, &Kr, &Fr, &u, &U, &S );
	return 0;
	
}


void Count_Number(int *np, int *ne){
	char buffer[256];
	FILE *fp;
	*np = 0;
	fp = fopen("point.txt", "r");
	while(fgets(buffer, 256, fp) != NULL){
		*np = *np + 1;
	}
	fclose(fp);
	*ne = 0;
	fp = fopen("tri_index.txt", "r");
	while(fgets(buffer, 256, fp) != NULL){
		*ne = *ne + 1;
	}
	fclose(fp);
//	*np = 11;
//	*ne = 12;
}


void Allocate_Memory(int np, int ne, char **Fix_x, char **Fix_y, int **tri \
		    , float **x, float **y, float **Fx, float **Fy, float **K \
		    , float **Kr, float **Fr, float **u, float **U, float **S){
	*Fix_x = ( char* )calloc( np , sizeof( char ) );
	*Fix_y = ( char* )calloc( np , sizeof( char ) );
	*tri = ( int* )malloc( 3 * ne * sizeof( int ) );
	*x = ( float* )malloc( np * sizeof( float ) );
	*y = ( float* )malloc( np * sizeof( float ) );
	*Fx = ( float* )calloc( np , sizeof( float ) );
	*Fy = ( float* )calloc( np , sizeof( float ) );
	*K = ( float* )calloc( 4 * np * np , sizeof( float ) );
	*Kr = ( float* )malloc( 4 * np * np * sizeof( float ) );
	*Fr = ( float* )malloc( 2 * np * sizeof( float ) );
	*u = ( float* )malloc( 2 * np * sizeof( float ) );
	*U = ( float* )malloc( 2 * np * sizeof( float ) );
	*S = ( float* )malloc( 2 * np * sizeof( float ) );

}

void Initial(int np, int ne, char *Fix_x, char *Fix_y, int *tri, float *x, float *y, float *Fx, float *Fy){
	int i;
	FILE *fp;
	fp = fopen( "point.txt" , "r" );
	for( i = 0 ; i < np ; ++i ){
		fscanf( fp, "%f %f\n", x + i , y + i );
		if( x[ i ] == 0.0 ){
			Fix_x[ i ] = 1;
			Fix_y[ i ] = 1;
		}
	}
	fclose(fp);
	fp = fopen( "tri_index.txt" , "r" );
	for( i = 0 ; i < ne ; ++i ){
		fscanf( fp, "%d %d %d\n", tri + 3 * i , tri + 1 + 3 * i , tri + 2 + 3 * i );
	}
	fclose(fp);
	Fy[ 3 ] = -12.5e3;

/*	x[0]=0;
       	x[1]=0.25;
       	x[2]=0.125; 
	x[3]=0; 
	x[4]=0.25; 
	x[5]=0.125; 
	x[6]=0; 
	x[7]=0.25; 
	x[8]=0.375; 
	x[9]=0.5; 
	x[10]=0.5;
	y[0]=0.5; 
	y[1]=0.5;
       	y[2]=0.375; 
	y[3]=0.25; 
	y[4]=0.25; 
	y[5]=0.125; 
	y[6]=0; 
	y[7]=0; 
	y[8]=0.125; 
	y[9]=0.25; 
	y[10]=0;
	Fy[4]=-12.5e3;
	Fix_x[0]=1;
	Fix_x[3]=1;
	Fix_x[6]=1;
	Fix_y[0]=1;
	Fix_y[3]=1;
	Fix_y[6]=1;
	//ele(1,:) = [1, 3, 2];
	tri[0]=0 ;tri[1]=2 ;tri[2]=1 ;
	//ele(2,:) = [1, 4, 3];
	tri[3]=0 ;tri[4]=3 ;tri[5]=2 ;
	//ele(3,:) = [3, 5, 2];
	tri[6]=2 ;tri[7]=4 ;tri[8]=1 ;
	//ele(4,:) = [3, 4, 5];
	tri[9]=2 ;tri[10]=3 ;tri[11]=4 ;
	//ele(5,:) = [4, 6, 5];
	tri[12]=3 ;tri[13]=5 ;tri[14]=4 ;
	//ele(6,:) = [4, 7, 6];
	tri[15]=3 ;tri[16]=6 ;tri[17]=5 ;
	//ele(7,:) = [5, 6, 8];
	tri[18]=4 ;tri[19]=5 ;tri[20]=7 ;
	//ele(8,:) = [6, 7, 8];
	tri[21]=5 ;tri[22]=6 ;tri[23]=7 ;
	//ele(9,:) = [5, 8 ,9];
	tri[24]=4 ;tri[25]=7 ;tri[26]=8 ;
	//ele(10,:) = [5, 9, 10];
	tri[27]=4 ;tri[28]=8 ;tri[29]=9 ;
	//ele(11,:) = [8, 11, 9];
	tri[30]=7 ;tri[31]=10 ;tri[32]=8 ;
	//ele(12,:) = [9, 11, 10];
	tri[33]=8 ;tri[34]=10 ;tri[35]=9 ;
*/
}

void Calculate_K(int np, int ne, int *tri, float *x, float *y, float *K){
	float A;
	float kl[36], D[9] , B[18], temp[ 18 ];
	int i,j,k,l;
	D[ 0 ] = 1  * E / ( 1 - NU * NU );
	D[ 1 ] = NU * E / ( 1 - NU * NU );
	D[ 2 ] = 0;

	D[ 3 ] = NU * E / ( 1 - NU * NU );
	D[ 4 ] = 1  * E / ( 1 - NU * NU );
	D[ 5 ] = 0;

	D[ 6 ] = 0;
	D[ 7 ] = 0;
	D[ 8 ] = 0.5 * ( 1 - NU ) * E / ( 1 - NU * NU );
	
	for( i = 0 ; i < ne ; ++i ){
		//calculate local k
		A = 0.5 * ( x[ tri[     i * 3] ] * ( y[ tri[ 1 + i * 3 ] ] - y[ tri[ 2 + i * 3 ] ] )\
			  + x[ tri[ 1 + i * 3] ] * ( y[ tri[ 2 + i * 3 ] ] - y[ tri[   + i * 3 ] ] )\
			  + x[ tri[ 2 + i * 3] ] * ( y[ tri[   + i * 3 ] ] - y[ tri[ 1 + i * 3 ] ] ) );
		B[ 0  ] = 0.5 / A * ( y[ tri[ 1 + i * 3 ] ] - y[ tri[ 2 + i * 3 ] ] );
		B[ 1  ] = 0;
		B[ 2  ] = 0.5 / A * ( y[ tri[ 2 + i * 3 ] ] - y[ tri[   + i * 3 ] ] );
		B[ 3  ] = 0;
		B[ 4  ] = 0.5 / A * ( y[ tri[   + i * 3 ] ] - y[ tri[ 1 + i * 3 ] ] );
		B[ 5  ] = 0;

		B[ 6  ] =  0;
		B[ 7  ] = -0.5 / A * ( x[ tri[ 1 + i * 3 ] ] - x[ tri[ 2 + i * 3 ] ] );
		B[ 8  ] =  0;
		B[ 9  ] = -0.5 / A * ( x[ tri[ 2 + i * 3 ] ] - x[ tri[   + i * 3 ] ] );
		B[ 10 ] =  0;
		B[ 11 ] = -0.5 / A * ( x[ tri[   + i * 3 ] ] - x[ tri[ 1 + i * 3 ] ] );

		B[ 12 ] = -0.5 / A * ( x[ tri[ 1 + i * 3 ] ] - x[ tri[ 2 + i * 3 ] ] );
		B[ 13 ] =  0.5 / A * ( y[ tri[ 1 + i * 3 ] ] - y[ tri[ 2 + i * 3 ] ] );
		B[ 14 ] = -0.5 / A * ( x[ tri[ 2 + i * 3 ] ] - x[ tri[   + i * 3 ] ] );
		B[ 15 ] =  0.5 / A * ( y[ tri[ 2 + i * 3 ] ] - y[ tri[   + i * 3 ] ] );
		B[ 16 ] = -0.5 / A * ( x[ tri[     i * 3 ] ] - x[ tri[ 1 + i * 3 ] ] );
		B[ 17 ] =  0.5 / A * ( y[ tri[     i * 3 ] ] - y[ tri[ 1 + i * 3 ] ] );
	
		for( j = 0 ; j < 3 ; ++j ){
			for( k = 0 ; k < 6 ; ++k ){
				temp[ k + j * 6 ] = 0;
				for( l = 0 ; l < 3 ; ++l ){
					temp[ k + j * 6 ] += D[ l + j * 3 ] * B[ k + l * 6 ];
				}
			}
		}

		for( j = 0 ; j < 6 ; ++j ){
			for( k = 0 ; k < 6 ; ++k ){
				kl[ k + j * 6 ] = 0;
				for( l = 0 ; l < 3 ; ++l ){
					kl[ k + j * 6 ] += B[ j + l * 6 ] * temp[ k + l * 6 ];
				}
			}
		}
		for( j = 0 ; j < 36 ; ++j ){
			kl[ j ] = kl[ j ] * A * THICK;
		}
		//Add into global k
		for( j = 0 ; j < 3 ; ++j ){
			for( k = 0 ; k < 3 ; ++k ){
				K[tri[k+i*3]*2  +(tri[j+i*3])*4*np     ]+=kl[k*2  +j*12  ];
				K[tri[k+i*3]*2+1+(tri[j+i*3])*4*np     ]+=kl[k*2+1+j*12  ];
				K[tri[k+i*3]*2  +(tri[j+i*3])*4*np+2*np]+=kl[k*2  +j*12+6];
				K[tri[k+i*3]*2+1+(tri[j+i*3])*4*np+2*np]+=kl[k*2+1+j*12+6];
			}
		}
	}
}

void Partition_K_F( int np, int *nr, char *Fix_x, char *Fix_y, float *Fx, float *Fy \
		   , float *Fr, float *K, float *Kr){
	int i , j , index = 0 , ix , iy;
	char *Fix;
	Fix = (char*)malloc( 2 * np * sizeof(char) );
	for( i = 0 ; i < np ; ++i ){
		Fix[2*i] = Fix_x[i];
		Fix[2*i+1] = Fix_y[i];
		if(Fix_x[i] == 0){
			Fr[index] = Fx[i];
			index++;
		}
		if(Fix_y[i] == 0 ){
			Fr[index] = Fy[i];
			index++;
		}
	}
	*nr = index;
	iy = 0;
	for( j = 0 ; j < 2*np ; ++j){
		ix = 0;
		if( Fix[ j ] == 0 ){
			for( i = 0 ; i < 2 * np ; ++i ){
				if( Fix[ i ] == 0 ){
					Kr[ ix + iy * index ] = K[ i + j * 2 * np ];
					ix++;
				}
			}
			iy++;
		}
	}
	free(Fix);
}

void Calc_Stress(){
	
}

void Save_Result(int np, int ne, int *tri, float *a, float *b){
	FILE *fp1, *fp2;
	int i,j;
	fp1 = fopen("data_x.txt", "w");
	fp2 = fopen("data_y.txt", "w");
	for(i = 0; i < ne; i++){
		fprintf(fp1,"%e %e %e\n", a[tri[3*i]], a[tri[1+3*i]], a[tri[2+3*i]]);	
		fprintf(fp2,"%e %e %e\n", b[tri[3*i]], b[tri[1+3*i]], b[tri[2+3*i]]);
	}
	fclose(fp1);
	fclose(fp2);
}


void Partition_U(int np, int nr, char *Fix_x , char *Fix_y , float *u, float *U){
	int il, ig, i;
	il=0;
	ig=0;
	for( i = 0 ; i < np ; i++ ){
		if(Fix_x[i] == 0){
			U[ig]=u[il];
			il++;
		}else{
			U[ig]=0.0;
		}
		ig++;
		if(Fix_y[i] == 0){
			U[ig]=u[il];
			il++;
		}else{
			U[ig]=0.0;
		}
		ig++;
	}
}

void Free_Memory( char **Fix_x, char **Fix_y, int **tri \
	 , float **x, float **y, float **Fx, float **Fy, float **K \
	 , float **Kr, float **Fr, float **u ,float **U, float **S ){
	free( *Fix_x );
	free( *Fix_y );
	free( *tri );
	free( *x );
	free( *y );
	free( *Fx );
	free( *Fy );
	free( *K );
	free( *Kr );
	free( *Fr );
	free( *u );
	free( *U );
	free( *S );
}

