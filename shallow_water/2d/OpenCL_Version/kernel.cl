#define G 9.81

__kernel void GPU_Calc(
		       int Nx,
		       int Ny,
		       float Zx,
		       float Zy,
		       __global float *A,
		       __global float *B,
		       __global float *C)
{
	float velx , vely , acc , F1 , F2 , F3 , Frx , Fry \
	    , FL1x , FL2x , FL3x , FR1x , FR2x , FR3x \
	    , FL1y , FL2y , FL3y , FR1y , FR2y , FR3y;
	int i = get_global_id( 0 );
	if( i < ( Nx + 2 ) * ( Ny + 2 ) ){
		velx = A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] / A[ i ];
		vely = A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] / A[ i ];
		acc    = sqrt( G * A[ i  ] );
		Frx = velx / acc;
		Fry = vely / acc;
		if (Frx >  1 ) Frx =  1;
		if (Frx < -1 ) Frx = -1;
		if (Fry >  1 ) Fry =  1;
		if (Fry < -1 ) Fry = -1;
		//x-dir
		F1 = A[ i ] * velx;
		F2 = A[ i ] * velx * velx + 0.5 * G * A[ i ] * A[ i ];
		F3 = A[ i ] * velx * vely;
		B[ i                                ]\
			=  0.5 * ( F1 * ( Frx + 1 )\
			+ A[ i                               ] * acc * ( 1 - Frx * Frx ) );
		B[ i  + ( Nx + 2 ) * ( Ny + 2 )     ]\
			=  0.5 * ( F2 * ( Frx + 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Frx * Frx ) );
		B[ i  + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]\
			=  0.5 * ( F3 * ( Frx + 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
		C[ i                                ]\
			= -0.5 * ( F1 * ( Frx - 1 )\
			+ A[ i                               ] * acc * ( 1 - Frx * Frx ) );
		C[ i + ( Nx + 2 ) * ( Ny + 2 )     ]\
			= -0.5 * ( F2 * ( Frx - 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Frx * Frx ) );
 		C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ]
			= -0.5 * ( F3 * ( Frx - 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Frx * Frx ) );
		//y-dir
		F1 = A[ i ] * vely;
		F2 = A[ i ] * vely * velx;
		F3 = A[ i ] * vely * vely + 0.5 * G * A[ i ] * A[ i ];
		B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			=  0.5 * ( F1 * ( Fry + 1 )\
			+ A[ i                               ] * acc * ( 1 - Fry * Fry ) );
		B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			=  0.5 * ( F2 * ( Fry + 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			=  0.5 * ( F3 * ( Fry + 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );
		C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			= -0.5 * ( F1 * ( Fry - 1 )\
			+ A[ i                               ] * acc * ( 1 - Fry * Fry ) );
		C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			= -0.5 * ( F2 * ( Fry - 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] * acc * ( 1 - Fry * Fry ) );
		C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			= -0.5 * ( F3 * ( Fry - 1 )\
			+ A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] * acc * ( 1 - Fry * Fry ) );	
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	
	//compute u
	if( ( i / ( Nx + 2 ) > 0 ) && ( i / ( Nx + 2 ) <= Nx ) && ( i % ( Nx + 2 ) > 0 ) && ( i % ( Nx + 2 ) <= Nx ) ){
			FL1x = B[ i - 1                                ]\
			     + C[ i                                    ];
			FR1x = B[ i                                    ]\
			     + C[ i + 1                                ];
			FL2x = B[ i - 1 + ( Nx + 2 ) * ( Ny + 2 )      ]\
			     + C[ i + ( Nx + 2 ) * ( Ny + 2 )          ];
			FR2x = B[ i + ( Nx + 2 ) * ( Ny + 2 )          ]\
			     + C[ i + 1 + ( Nx + 2 ) * ( Ny + 2 )      ];
			FL3x = B[ i - 1 + ( Nx + 2 ) * ( Ny +  2 ) * 2 ]\
			     + C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2      ];
			FR3x = B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2      ]\
			     + C[ i + 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2  ];

			FL1y = B[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ]\
			     + C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3              ];
			FR1y = B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 3              ]\
			     + C[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 3 ];
			FL2y = B[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ]\
			     + C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4              ];
			FR2y = B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 4              ]\
			     + C[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 4 ];
			FL3y = B[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ]\
			     + C[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5              ];
			FR3y = B[ i + ( Nx + 2 ) * ( Ny + 2 ) * 5              ]\
			     + C[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 5 ];

			A[ i                              ]\
				= A[ i                              ]\
				- Zx * ( FR1x - FL1x ) - Zy * ( FR1y - FL1y );
			A[ i + ( Nx + 2 ) * ( Ny + 2)     ]\
				= A[ i + ( Nx + 2 ) * ( Ny + 2)     ]\
				- Zx * ( FR2x - FL2x ) - Zy * ( FR2y - FL2y );
			A[ i + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
				= A[ i + ( Nx + 2 ) * ( Ny + 2) * 2 ]\
				- Zx * ( FR3x - FL3x ) - Zy * ( FR3y - FL3y );	
	}
	barrier( CLK_LOCAL_MEM_FENCE );

	//set reflective boundary
	if( i % ( Nx + 2 ) == 0 ){
		A[ i                               ] \
		=   A[ i + 1                               ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		= - A[ i + 1 + ( Nx + 2 ) * ( Ny + 2 )     ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		=   A[ i + 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i % ( Nx + 2 ) == Nx + 1 ){
		A[ i                               ] \
		=   A[ i - 1                               ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		= - A[ i - 1 + ( Nx + 2 ) * ( Ny + 2 )     ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		=   A[ i - 1 + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i / ( Nx + 2 ) == 0 ){
		A[ i                               ] \
		=   A[ i + ( Nx + 2 )                               ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		=   A[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 )     ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		= - A[ i + ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}
	if( i / ( Nx + 2 ) == Ny + 1 ){
		A[ i                               ] \
		=   A[ i - ( Nx + 2 )                               ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 )     ] \
		=   A[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 )     ];
		A[ i + ( Nx + 2 ) * ( Ny + 2 ) * 2 ] \
		= - A[ i - ( Nx + 2 ) + ( Nx + 2 ) * ( Ny + 2 ) * 2 ];
	}

}
