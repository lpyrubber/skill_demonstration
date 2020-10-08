#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>

#define Lx 100.0
#define Ly 100.0
#define Nx 200
#define Ny 200
#define	DX (Lx/Nx)
#define DY (Ly/Ny)
#define DT (0.01*DX)
#define Zx (DT/DX)
#define Zy (DT/DY)
#define NO_STEP 1
#define G 9.81

float *u, *u_new;
float *fm, *fp;
FILE *pFile;

void Allocate_Memory();
void Initial();
void Compute_Flux();
void Compute_U();
void Update_U();
void Save_Result();
void Free_Memory();

int main(){
	static int k=0;
        char temp[20];
        int i;
        sprintf(temp, "data_%d.txt",k);
        pFile = fopen(temp, "w");

	Allocate_Memory();
	Initial();
	for( i = 0 ; i < NO_STEP ; ++i ){
		Compute_Flux();
		Compute_U();
		Update_U();
	}
	Save_Result();
	Free_Memory();
	fclose(pFile);
	return 0;
}

void Allocate_Memory(){
	size_t size = 3 * ( Nx + 2 ) * ( Ny + 2 ) * sizeof( float );
	u     = ( float* )malloc( size );
	u_new = ( float* )malloc( size );
	fm    = ( float* )malloc( 2 * size );
	fp    = ( float* )malloc( 2 * size );
}

void Initial(){
	int i , j;
	for( j = 0 ; j < Ny + 2 ; ++j ){
		for( i = 0 ; i < Nx + 2 ; ++i ){
//			if( ( ( i - 1 ) < 0.5 * Nx ) && ( ( j - 1 ) < 0.5 * Ny ) ){
			if( (  i - 1 ) < 0.5 * Nx ){
				u[ i + ( Nx + 2 ) * j                              ] = 10.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)     ] = 0.0;
				u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2) * 2 ] = 0.0;
			}else{
				u[ i + ( Nx + 2 ) * j                              ] = 1.0;
				u[ i + ( Nx + 2 ) * j + ( Ny + 2 ) * ( Ny + 2)     ] = 0.0;
				u[ i + ( Nx + 2 ) * j + ( Ny + 2 ) * ( Ny + 2) * 2 ] = 0.0;
			}
		}
	}
//	printf("finish Initial\n");
}

void Compute_Flux(){
	int i , j;
	float velx , vely , a;
	float F1 , F2 , F3;
	float Frx, Fry;
	for( j = 0 ; j < Ny + 2 ; ++j ){
		for( i = 0 ; i < Nx + 2 ; ++i ){
	                velx = u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
			       / u[ i + ( Nx + 2 ) * j ];
			vely = u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			       / u[ i + ( Nx + 2 ) * j ];
	                a    = sqrt( G * u[ i + ( Nx + 2 ) * j ] );
	                Frx = velx / a;
			Fry = vely / a;
	                if (Frx >  1 ) Frx =  1;
	                if (Frx < -1 ) Frx = -1;
	                if (Fry >  1 ) Fry =  1;
	                if (Fry < -1 ) Fry = -1;
			//x-dir
	                F1 = u[ i + ( Nx + 2 ) * j ] * velx;
	                F2 = u[ i + ( Nx + 2 ) * j ] * velx * velx\
			   + 0.5 * G * u[ i + ( Nx + 2 ) * j ] * u[ i + ( Nx + 2 ) * j ];
			F3 = u[ i + ( Nx + 2 ) * j ] * velx * vely;
	                fp[ i + ( Nx + 2 ) * j                               ]\
			       	=  0.5 * ( F1 * ( Frx + 1 )\
				+ u[ i + ( Nx + 2 ) * j                               ] * a * ( 1 - Frx * Frx ) );
	                fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)      ]\
			   	=  0.5 * ( F2 * ( Frx + 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] * a * ( 1 - Frx * Frx ) );
	                fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			       	=  0.5 * ( F2 * ( Frx + 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * a * ( 1 - Frx * Frx ) );
	                fm[ i + ( Nx + 2 ) * j                               ]\
			       	= -0.5 * ( F1 * ( Frx - 1 )\
				+ u[ i + ( Nx + 2 ) * j                               ] * a * ( 1 - Frx * Frx ) );
			fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
			       	= -0.5 * ( F2 * ( Frx - 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] * a * ( 1 - Frx * Frx ) );
                	fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]
			       	= -0.5 * ( F1 * ( Frx - 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * a * ( 1 - Frx * Frx ) );
			//y-dir
			F1 = u[ i + ( Nx + 2 ) * j ] * vely;
	                F2 = u[ i + ( Nx + 2 ) * j ] * vely * velx;
			F3 = u[ i + ( Nx + 2 ) * j ] * vely * vely\
			   + 0.5 * G * u[ i + ( Nx + 2 ) * j ] * u[ i + ( Nx + 2 ) * j ];
        	        fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			       	=  0.5 * ( F1 * ( Fry + 1 )\
				+ u[ i + ( Nx + 2 ) * j                               ] * a * ( 1 - Fry * Fry ) );
	                fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
				=  0.5 * ( F2 * ( Fry + 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] * a * ( 1 - Fry * Fry ) );
	                fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
				=  0.5 * ( F2 * ( Fry + 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * a * ( 1 - Fry * Fry ) );
	                fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			       	= -0.5 * ( F1 * ( Fry - 1 )\
				+ u[ i + ( Nx + 2 ) * j                               ] * a * ( 1 - Fry * Fry ) );
	                fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			       	= -0.5 * ( F2 * ( Fry - 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ] * a * ( 1 - Fry * Fry ) );
	                fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			       	= -0.5 * ( F1 * ( Fry - 1 )\
				+ u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * a * ( 1 - Fry * Fry ) );
		}
        }
//	printf("finish compute flux\n");
}

void Compute_U(){
	float FL1x , FL2x , FL3x ,  FR1x , FR2x , FR3x;
	float FL1y , FL2y , FL3y ,  FR1y , FR2y , FR3y;
	int i , j;
	for( j = 1 ; j < Ny + 1 ; ++j ){
        	for( i = 1 ; i < Nx + 1 ; i++ ){
	                FL1x = fp[ i - 1 + ( Nx + 2 ) * j                               ]\
			     + fm[ i + ( Nx + 2 ) * j                                   ];
	                FR1x = fp[ i + ( Nx + 2 ) * j                                   ]\
			     + fm[ i + 1 + ( Nx + 2 ) * j                               ];
	                FL2x = fp[ i - 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
			     + fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )         ];
	                FR2x = fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )         ]\
			     + fm[ i + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ];
			FL3x = fp[ i - 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			     + fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2     ];
			FR3x = fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2     ]\
			     + fm[ i + 1 + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];

	                FL1y = fp[ i + ( Nx + 2 ) * ( j - 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			     + fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 3         ];
	                FR1y = fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 3         ]\
			     + fm[ i + ( Nx + 2 ) * ( j + 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ];
	                FL2y = fp[ i + ( Nx + 2 ) * ( j - 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			     + fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 4        ];
	                FR2y = fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 4        ]\
			     + fm[ i + ( Nx + 2 ) * (j + 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ];
			FL3y = fp[ i + ( Nx + 2 ) * (j - 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			     + fm[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 5        ];
			FR3y = fp[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 5        ]\
			     + fm[ i + ( Nx + 2 ) * (j + 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ];

		       	u_new[ i + ( Nx + 2 ) * j                              ]\
				= u[ i + ( Nx + 2 ) * j                              ]\
			          - Zx * ( FR1x - FL1x ) + Zy * ( FR1y - FL1y );
	                u_new[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)     ]\
			       	= u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2)     ]\
				  - Zx * ( FR2x - FL2x ) + Zy * ( FR2y - FL2y );
	                u_new[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
			       	= u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
				  - Zx * ( FR3x - FL3x ) + Zy * ( FR3y - FL3y );
		}
        }
//	printf("\tfinish phase 1\n");
	for( i = 0 ; i < Nx + 2 ; ++i){
		u_new[ i                                                         ]\
		       	=   u_new[ i + Nx + 2                                       ];
		u_new[ i + ( Nx + 2 ) * ( Ny + 1 )                               ]\
			=   u_new[ i + ( Nx + 2 ) * Ny                              ];
		u_new[ i + ( Nx + 2 ) * ( Ny + 2 )                               ]\
			= - u_new[ i + Nx + 2 + ( Nx + 2 ) * ( Ny + 2 )             ];
		u_new[ i + ( Nx + 2 ) * ( Ny + 1 ) + ( Nx + 2 ) * ( Ny + 2 )     ]\
			= - u_new[ i + ( Nx + 2 ) * Ny + ( Nx + 2 ) * ( Ny + 2 )    ];
		u_new[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2                           ]\
			= - u_new[ i + Nx + 2 + ( Nx + 2 ) * ( Ny + 2 ) *2          ];
		u_new[ i + ( Nx + 2 ) * ( Ny + 1 ) + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			= - u_new[ i + ( Nx + 2 ) * Ny + ( Nx + 2 ) * ( Ny + 2 ) *2 ];
	}
	for( j = 0 ; j < Ny + 2 ; ++j){
		u_new[ j * Nx                                        ]\
			= u_new[ 1 + j * Nx                                ];
		u_new[ Nx + 1 + j * Nx                               ]\
			= u_new[ Nx + j * Nx                               ]; 
		u_new[ j * Nx + ( Nx + 2 ) * ( Ny + 2 )              ]\
			= u_new[ 1 + j * Nx + ( Nx + 2 ) * ( Ny + 2 )      ];
		u_new[ Nx + 1 + j * Nx + ( Nx + 2 ) * ( Ny + 2 )     ]\
			= u_new[ Nx + j * Nx + ( Nx + 2 ) * ( Ny + 2 )     ]; 
		u_new[ j * Nx + ( Nx + 2 ) * ( Ny + 2 ) * 2          ]\
			= u_new[ 1 + j * Nx + ( Nx + 2 ) * ( Ny + 2 ) * 2  ];
		u_new[ Nx + 1 + j * Nx + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			= u_new[ Nx + j * Nx + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]; 
	}
//	printf("finish compute u\n");
}

void Update_U(){
        int i, j;
	for(j = 0 ; j < Ny + 2 ; ++j){
	        for(i = 0 ; i < Nx +2 ; i++){
        	        u[ i + ( Nx + 2 ) * j                               ]\
				= u_new[ i + ( Nx + 2 ) * j                               ];
        	        u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
			       	= u_new[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ];
        	        u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			       	= u_new[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
		}
	}
//	printf("finish update u\n");
}

void Save_Result(){
        int i, j;
        if (pFile == NULL){
                printf("File open failed\n");
        }else{
		for( j = 1 ; j < Ny + 1 ; ++j){
 	        	for( i = 1 ; i < Nx + 1 ; ++i){
	        		fprintf(pFile, "%e %e %e %e %e\n"\
					, i * DX\
					, j * DY\
					, u[ i + ( Nx + 2 ) * j                               ]\
					, u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 )     ]\
					, u[ i + ( Nx + 2 ) * j + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]);
        	        }
		}
        }
//	printf("finish save result\n");
}

void Free_Memory(){
	free(fm);
	free(fp);
	free(u);
	free(u_new);
}
