#include "main.h"

void Allocate_Memory( int id , int max , float **a , float **b , float **c ){
	size_t size;
	int averow , extra , rows;
	if( id == MASTER){
		size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
		*a = ( float* )malloc( size );
		*b = ( float* )malloc( sizeof(float) );
		*c = ( float* )malloc( sizeof(float) );
	}else{
		averow = Ny / ( max - 1 );
		extra  = Ny % ( max - 1 );
		// +2 for ghost blocks
		rows = ( id <= extra ) ? ( averow + 3 ) * ( Nx + 2 ) : ( averow + 2 ) * (Nx + 2 );
		size = 3 * ( Nx + 2 ) * rows * sizeof( float );
		*a = ( float* )malloc( size );
		*b = ( float* )malloc( 2 * size );
		*c = ( float* )malloc( 2 * size );		
	}
}

void Compute( int id , int ele , float *a , float *b , float *c ){
	if( id != MASTER ){
		int i , j;
		float velx , vely , acc;
		float F1 , F2 , F3;
		float Frx, Fry;
		float FL1x , FL2x , FL3x ,  FR1x , FR2x , FR3x;
		float FL1y , FL2y , FL3y ,  FR1y , FR2y , FR3y;
		//compute flux
		for( j = 0 ; j < ele + 2 ; ++j ){
			for( i = 0 ; i < Nx + 2 ; ++i ){
		                velx = a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] / a[ i + ( Nx + 2 ) * j ];
				vely = a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] / a[ i + ( Nx + 2 ) * j ];
		                acc  = sqrt( G * a[ i + ( Nx + 2 ) * j ] );
		                Frx = velx / acc;
				Fry = vely / acc;
		                if (Frx >  1 ) Frx =  1;
		                if (Frx < -1 ) Frx = -1;
		                if (Fry >  1 ) Fry =  1;
		                if (Fry < -1 ) Fry = -1;
				//x-dir
		                F1 = a[ i + ( Nx + 2 ) * j ] * velx;
		                F2 = a[ i + ( Nx + 2 ) * j ] * velx * velx \
				   + 0.5 * G * a[ i + ( Nx + 2 ) * j ] * a[ i + ( Nx + 2 ) * j ];
				F3 = a[ i + ( Nx + 2 ) * j ] * velx * vely;
		                b[ i + ( Nx + 2 ) * j                                ]\
				       	=  0.5 * ( F1 * ( Frx + 1 )\
					+ a[ i + ( Nx + 2 ) * j                                ] * acc * ( 1 - Frx * Frx ) );
		                b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ]\
				   	=  0.5 * ( F2 * ( Frx + 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] * acc * ( 1 - Frx * Frx ) );
		                b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ]\
				       	=  0.5 * ( F3 * ( Frx + 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
		                c[ i + ( Nx + 2 ) * j                                ]\
				       	= -0.5 * ( F1 * ( Frx - 1 )\
					+ a[ i + ( Nx + 2 ) * j                                ] * acc * ( 1 - Frx * Frx ) );
				c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ]\
				       	= -0.5 * ( F2 * ( Frx - 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] * acc * ( 1 - Frx * Frx ) );
      		         	c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ]
				       	= -0.5 * ( F3 * ( Frx - 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
				//y-dir
				F1 = a[ i + ( Nx + 2 ) * j ] * vely;
		                F2 = a[ i + ( Nx + 2 ) * j ] * vely * velx;
				F3 = a[ i + ( Nx + 2 ) * j ] * vely * vely\
				   + 0.5 * G * a[ i + ( Nx + 2 ) * j ] * a[ i + ( Nx + 2 ) * j ];
        		        b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 3 ]\
				       	=  0.5 * ( F1 * ( Fry + 1 )\
					+ a[ i + ( Nx + 2 ) * j                                ] * acc * ( 1 - Fry * Fry ) );
		                b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 4 ]\
					=  0.5 * ( F2 * ( Fry + 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		                b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 5 ]\
					=  0.5 * ( F3 * ( Fry + 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );
		                c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 3 ]\
				       	= -0.5 * ( F1 * ( Fry - 1 )\
					+ a[ i + ( Nx + 2 ) * j                                ] * acc * ( 1 - Fry * Fry ) );
		                c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 4 ]\
				       	= -0.5 * ( F2 * ( Fry - 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		                c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 5 ]\
				       	= -0.5 * ( F3 * ( Fry - 1 )\
					+ a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );
					
			}
		}
		//compute u
		for( j = 1 ; j < ele + 1 ; ++j ){
	        	for( i = 1 ; i < Nx + 1 ; i++ ){
		                FL1x = b[ i - 1 + ( Nx + 2 ) * j                                ]\
				     + c[ i + ( Nx + 2 ) * j                                    ];
		                FR1x = b[ i + ( Nx + 2 ) * j                                    ]\
				     + c[ i + 1 + ( Nx + 2 ) * j                                ];
		                FL2x = b[ i - 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ]\
				     + c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )         ];
		                FR2x = b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )         ]\
				     + c[ i + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 )     ];
				FL3x = b[ i - 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ]\
				     + c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2     ];
				FR3x = b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2     ]\
				     + c[ i + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 2 ];
	
		                FL1y = b[ i + ( Nx + 2 ) * ( j - 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 3 ]\
				     + c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 3         ];
		                FR1y = b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 3         ]\
				     + c[ i + ( Nx + 2 ) * ( j + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 3 ];
		                FL2y = b[ i + ( Nx + 2 ) * ( j - 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 4 ]\
				     + c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 4        ];
		                FR2y = b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 4        ]\
				     + c[ i + ( Nx + 2 ) * (j + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 4 ];
				FL3y = b[ i + ( Nx + 2 ) * (j - 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 5 ]\
				     + c[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 5        ];
				FR3y = b[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2 ) * 5        ]\
				     + c[ i + ( Nx + 2 ) * (j + 1 ) + ( Nx + 2 ) * ( ele + 2 ) * 5 ];
//				if(j==1)printf("%d %f %f %f %f %f %f\n", i, FL1x, FR1x, FL2x, FR2x, FL3x, FR3x);
			       	a[ i + ( Nx + 2 ) * j                    ]\
					= a[ i + ( Nx + 2 ) * j                              ]\
				          - Zx * ( FR1x - FL1x ) - Zy * ( FR1y - FL1y );
		                a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2)     ]\
				       	= a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2)     ]\
					  - Zx * ( FR2x - FL2x ) - Zy * ( FR2y - FL2y );
		                a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2) * 2 ]\
				       	= a[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( ele + 2) * 2 ]\
					  - Zx * ( FR3x - FL3x ) - Zy * ( FR3y - FL3y );
			}
	        }
	}
	
}
void Free_Memory( float **a , float **b , float **c ){
	free(*a);
	free(*b);
	free(*c);
}

