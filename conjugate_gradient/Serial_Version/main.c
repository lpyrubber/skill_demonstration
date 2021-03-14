#include "main.h"

int main(){
	float *A,  *b,  *x;
	Allocate_Memory( &A , &b , &x );
	Initial( A , b , x );
	Conjugate_Gradient( A , b , x );
	Save_Result( x );
	Free_Memory( &A , &b , &x );
}

void Allocate_Memory( float **A , float **b, float **x ){
	size_t size;
	size = N * N * sizeof( float );
	*A = ( float* )malloc( size );
	size = N * sizeof( float );
	*b = ( float* )malloc( size );
	*x = ( float* )malloc( size );
}

void Initial( float *A , float *b , float *x ){
	int i , j;
	for( j = 0 ; j < N ; ++j ){
		x[ j ] = 0;
		b[ j ] = ( j == 0 || j == N - 1 ) ? 0 : Q/K*dx*dx ;
		for( i = 0 ; i < N ; ++i ){
			if(i == j){
				A[ i + N * j ] = 2;
			}else if( j == i - 1 && j>0 && j<N-1 ){
				A[ i + N * j ] = -1;
			}else if( j == i + 1 && j>0 && j<N-1 ){
				A[ i + N * j ] = -1;
			}else{
				A[ i + N * j ] = 0;
			}
		}
	}
}

void Conjugate_Gradient( float *A , float *b, float *x ){
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
			r[ j ] -= A[ i + N * j ] * x[ i ];
		}
		p[ j ] = r[ j ];
		rs_old += r[ j ] * r[ j ];
	}
	k = 0;
	flag = 1;
	while( flag && k < 2){
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
		if( sqrt( rs_new ) < 1e-10 ){
			flag = 0;
			printf("k=%d\n", k);
		}
		for( i = 0 ; i < N ; ++i ){
			p[ i ] = r[ i ] + ( rs_new / rs_old ) * p[ i ];
		}
		rs_old = rs_new;
		printf("rs_old = %e\n", rs_old );
		k++;
	}
	free(p);
	free(s);
	free(r);
}

void Save_Result( float *x ){
	FILE *fp;
	int i;
	fp = fopen( "data.txt", "w" );
	for( i = 0 ; i < N ; ++i ){
		fprintf( fp , "%f\n", x[ i ] );
	}
}

void Free_Memory( float **A , float **b , float **x ){
	free(*A);
	free(*b);
	free(*x);
}
