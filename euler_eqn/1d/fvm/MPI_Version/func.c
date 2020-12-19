#include "main.h"

void Allocate_Memory( int id , int max , float **a , float **b , float **c , float **d , float **e , float **f){
	//define variable in function
	int i , averow , extra , rows ;
	averow = N / max;
	extra  = N % max;
	if( id == MASTER ){
		*a = ( float* )malloc( 3 * (N+1) * sizeof( float ) );
		*b = ( float* )malloc( (N+1) * sizeof( float ) );
		*c = ( float* )malloc( (N+1) * sizeof( float ) );
		*d = ( float* )malloc( (N+1) * sizeof( float ) );
		*e = ( float* )malloc( 3 * ( averow + 3 ) * sizeof( float ) );
		*f = ( float* )malloc( 3 * ( averow + 3 ) * sizeof( float ) );
	}else{
		rows = ( id < extra ) ?  averow + 3 : averow + 2 ; // +2 for ghost blocks
		*a = ( float* )malloc( 3 * rows * sizeof( float ) );
		*b = ( float* )malloc( rows * sizeof( float ) );
		*c = ( float* )malloc( rows * sizeof( float ) );
		*d = ( float* )malloc( rows * sizeof( float ) );
		*e = ( float* )malloc( 3 * rows * sizeof( float ) );	
		*f = ( float* )malloc( 3 * rows * sizeof( float ) );	
	}
}

void Compute( int id , int ele , float *a , float *b , float *c , float *d , float *e , float *f ){
	//create variable in this function
	int i, len;
	float acc , F1 , F2 , F3 , Fr , FL1 , FL2 , FL3 , FR1 , FR2 , FR3;
	
	//compute flux
	len = (id == 0)? N+1 : ele + 2;
	for(i =0; i < ele + 2 ; i++ ){
		acc = sqrt(gamma*R*d[i]);
		Fr = c[i]/acc;
		if (Fr >1 ) Fr =1;
		if (Fr < -1) Fr = -1;
		F1 = a[i]*c[i];
		F2 = a[i]*c[i]*c[i] + R*b[i]*d[i];
		F3 = c[i]*(a[i+len*2] + R*b[i]*d[i]);
		f[i                   ] = 0.5*(F1*(Fr+1)+a[i]*acc*(1-Fr*Fr));
		f[i + ele + 2         ] = 0.5*(F2*(Fr+1)+a[i+len]*acc*(1-Fr*Fr));
		f[i + ( ele + 2 ) * 2 ] = 0.5*(F3*(Fr+1)+a[i+len*2]*acc*(1-Fr*Fr));
		e[i                   ] = -0.5*(F1*(Fr-1)+a[i]*acc*(1-Fr*Fr));
		e[i + ele + 2         ] = -0.5*(F2*(Fr-1)+a[i+len]*acc*(1-Fr*Fr));
		e[i + ( ele + 2 ) * 2 ] = -0.5*(F3*(Fr-1)+a[i+len*2]*acc*(1-Fr*Fr));
	}
	//compute u
	for(i=1; i< ele + 1; i++){
		FL1 = f[ i - 1 ]+e[ i ];
		FR1 = f[ i ]+e[ i + 1 ];
		FL2 = f[ i - 1 + ele + 2]+e[ i + ele + 2 ];
		FR2 = f[ i + ele + 2]+e[ i + 1 + ele + 2 ];
		FL3 = f[ i - 1 + ( ele + 2 ) * 2 ]+e[ i + ( ele + 2 ) * 2 ];
		FR3 = f[ i + ( ele + 2 ) * 2 ]+e[ i + 1 + ( ele + 2 ) * 2 ];
		a[ i           ] = a[ i           ] - Z*(FR1-FL1);
		a[ i + len     ] = a[ i + len     ] - Z*(FR2-FL2);
		a[ i + len * 2 ] = a[ i + len * 2 ] - Z*(FR3-FL3);
	}
	for(i=0; i< ele + 2; ++i){
		b[i] = a[i];
		c[i] = a[i+len]/a[i];
		d[i] = (a[i+len*2]/a[i]-c[i]*c[i]*0.5)/CV;
	}
}

void Free_Memory( float **a , float **b , float **c , float **d , float **e , float **f ){
	free( *a );
	free( *b );
	free( *c );
	free( *d );
	free( *e );
	free( *f );	
}
