#include "main.h"

void Allocate_Memory( int id , int max , float **a , float **b , float **c , float **d , float **e , float **f , float **g , float **h , float **i){
	//define variable in function
	int  averow , extra , rows ;
	averow = N / max;
	extra  = N % max;
	if( id == MASTER ){
		*a = ( float* )malloc( 3 * (N+1) * sizeof( float ) );
		*b = ( float* )malloc( (N+1) * sizeof( float ) );
		*c = ( float* )malloc( (N+1) * sizeof( float ) );
		*d = ( float* )malloc( (N+1) * sizeof( float ) );
		*e = ( float* )malloc( 3 * ( averow + 3 ) * sizeof( float ) );
		*f = ( float* )malloc( ( averow + 3 ) * sizeof( float ) );
		*g = ( float* )malloc( ( averow + 3 ) * sizeof( float ) );
		*h = ( float* )malloc( ( averow + 3 ) * sizeof( float ) );
		*i = ( float* )malloc( ( averow + 3 ) * sizeof( float ) );
	}else{
		rows = ( id < extra ) ?  averow + 3 : averow + 2 ; // +2 for ghost blocks
		*a = ( float* )malloc( 3 * rows * sizeof( float ) );
		*b = ( float* )malloc( rows * sizeof( float ) );
		*c = ( float* )malloc( rows * sizeof( float ) );
		*d = ( float* )malloc( rows * sizeof( float ) );
		*e = ( float* )malloc( 3 * rows * sizeof( float ) );	
		*f = ( float* )malloc( rows * sizeof( float ) );	
		*g = ( float* )malloc( rows * sizeof( float ) );	
		*h = ( float* )malloc( rows * sizeof( float ) );	
		*i = ( float* )malloc( rows * sizeof( float ) );	
	}
}

void Compute( int id , int ele , float *a , float *b , float *c , float *d , float *e , float *f , float *g , float *h , float *i ){
	//create variable in this function
	int i1, len;
	float acc, M, Mt, P;
	float FL[3], FR[3];
	//compute flux
	len = (id == 0)? N+1 : ele + 2;
	for(i1 =0; i1 < ele + 2 ; i1++ ){
		acc=sqrt(gamma*R*d[i1]);
		M=c[i1]/acc;
		P=b[i1]*R*d[i1];
		Mt=fabs(M);
		if(Mt<=1){
			f[i1]= 0.25*P*(M+1)*(M+1)*(2-M);
			g[i1]= 0.25*P*(M-1)*(M-1)*(2+M);
			h[i1]= 0.25*(M+1)*(M+1);
			i[i1]=-0.25*(M-1)*(M-1);
		}else{
			f[i1]= 0.5*P*(M+Mt)/M;
			g[i1]= 0.5*P*(M-Mt)/M;
			h[i1]= 0.5*(M+Mt);
			i[i1]=-0.5*(M-Mt);
		}
		e[i1]=b[i1]*acc;
		e[i1+ele+2]=b[i1]*acc*c[i1];
		e[i1+(ele+2)*2]=b[i1]*acc*(CP*d[i1]+0.5*c[i1]*c[i1]);
	}
	//compute u
	for(i1=1; i1 < ele + 1; i1++){
		Mt=h[i1-1]+i[i1];
		FL[0]=Mt*(e[i1-1          ]+e[i1          ])*0.5-fabs(Mt)*(e[i1          ]-e[i1-1          ])*0.5;
		FL[1]=Mt*(e[i1-1+ele+2    ]+e[i1+ele+2    ])*0.5-fabs(Mt)*(e[i1+ele+2    ]-e[i1-1+ele+2    ])*0.5+(f[i1-1]+g[i1]);
		FL[2]=Mt*(e[i1-1+(ele+2)*2]+e[i1+(ele+2)*2])*0.5-fabs(Mt)*(e[i1+(ele+2)*2]-e[i1-1+(ele+2)*2])*0.5;
		Mt=h[i1]+i[i1+1];
		FR[0]=Mt*(e[i1          ]+e[i1+1          ])*0.5-fabs(Mt)*(e[i1+1          ]-e[i1          ])*0.5;
		FR[1]=Mt*(e[i1+ele+2    ]+e[i1+1+ele+2    ])*0.5-fabs(Mt)*(e[i1+1+ele+2    ]-e[i1+ele+2    ])*0.5+(f[i1]+g[i1+1]);
		FR[2]=Mt*(e[i1+(ele+2)*2]+e[i1+1+(ele+2)*2])*0.5-fabs(Mt)*(e[i1+1+(ele+2)*2]-e[i1+(ele+2)*2])*0.5;
		a[i1      ]-=Z*(FR[0]-FL[0]);
		a[i1+len  ]-=Z*(FR[1]-FL[1]);
		a[i1+len*2]-=Z*(FR[2]-FL[2]);
	}
	for(i1=0; i1< ele + 2; ++i1){
		b[i1] = a[i1];
		c[i1] = a[i1+len]/a[i1];
		d[i1] = (a[i1+len*2]/a[i1]-c[i1]*c[i1]*0.5)/CV;
	}
}

void Free_Memory( float **a , float **b , float **c , float **d , float **e , float **f , float **g , float **h , float **i ){
	free( *a );
	free( *b );
	free( *c );
	free( *d );
	free( *e );
	free( *f );	
	free( *g );	
	free( *h );	
	free( *i );	

}
