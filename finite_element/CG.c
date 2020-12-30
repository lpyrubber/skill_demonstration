#include "main.h"

void Conjugate_Gradient( float *A , float *b, float *x , int N){
	size_t size;
	float *p , *r, *s;
	float alpha, beta, rs_old, rs_new, temp;
	int i, j, k;
	char flag;
	size = N * sizeof( float );
	p = ( float* )malloc( size );
	r = ( float* )malloc( size );
	s = ( float* )malloc( size );
	
	rs_old = 0;
	for( j = 0 ; j < N ; ++j ){
		r[ j ] = b[ j ];
		for( i = 0 ; i < N ; ++i ){
			r[ j ] -= A[ i + N * j ] * x[ j ];
		}
		p[ j ] = r[ j ];
		rs_old += r[ j ] * r[ j ];
	}
	k = 0;
	flag = 1;
	while( flag && k < 1e5){
		temp = 0;
		for( j = 0 ; j < N ; ++j ){
			s[ j ] = 0;
			for( i = 0 ; i < N ; ++i ){
				s[ j ] += A[ i + N * j ] * p[ i ];
			}
			temp += s[ j ] * p[ j ];
		}
		alpha = rs_old / temp;
		rs_new = 0;
		for( i = 0 ; i < N ; ++i ){
			x[ i ] += alpha * p[ i ];
			r[ i ] -= alpha * s[ i ];
			rs_new += r[ i ] * r[ i ];
		}
		if( rs_new  < 1e-6 ){
			flag = 0;
			printf("k=%d\n", k);
		}
		for( i = 0 ; i < N ; ++i ){
			p[ i ] = r[ i ] + ( rs_new / rs_old ) * p[ i ];
		}
		rs_old = rs_new;
		k++;
	}
	free(p);
	free(s);
	free(r);
}
