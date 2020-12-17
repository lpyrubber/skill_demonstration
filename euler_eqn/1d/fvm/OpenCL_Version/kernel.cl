#define R 1.0
#define gamma 1.4
#define CV (R/(gamma-1))

__kernel void GPU_Calc(
		       int N,
		       float Z,
		       __global float *A,
		       __global float *B,
		       __global float *C,
		       __global float *D,
		       __global float *E,
		       __global float *F)
{
	float a, F1, F2, F3, Fr;
	float FL1, FL2, FL3, FR1, FR2, FR3;
	int i = get_global_id( 0 );
	if( i < N ){
		a = sqrt( gamma * R * D[ i ] );
		Fr = C[ i ] / a;
		if ( Fr > 1  ) Fr =  1;
		if ( Fr < -1 ) Fr = -1;
		F1 = A[ i ] * C[ i ];
		F2 = A[ i ] * C[ i ] * C[ i ] +  B[ i ] * R * D[ i ];
		F3 = C[ i ] * ( A[ i + N * 2 ] + B[ i ] * R * D[ i ] );
		F[ i         ] =  0.5*(F1*(Fr+1)+A[ i         ]*a*(1-Fr*Fr));
		F[ i + N     ] =  0.5*(F2*(Fr+1)+A[ i + N     ]*a*(1-Fr*Fr));
		F[ i + N * 2 ] =  0.5*(F3*(Fr+1)+A[ i + N * 2 ]*a*(1-Fr*Fr));
		E[ i         ] = -0.5*(F1*(Fr-1)+A[ i         ]*a*(1-Fr*Fr));
		E[ i + N     ] = -0.5*(F2*(Fr-1)+A[ i + N     ]*a*(1-Fr*Fr));
		E[ i + N * 2 ] = -0.5*(F3*(Fr-1)+A[ i + N * 2 ]*a*(1-Fr*Fr));
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if( i > 0 && i < N - 1 ){
		FL1 = F[ i - 1         ]+E[ i             ];
		FR1 = F[ i             ]+E[ i + 1         ];
		FL2 = F[ i - 1 + N     ]+E[ i     + N     ];
		FR2 = F[ i     + N     ]+E[ i + 1 + N     ];
		FL3 = F[ i - 1 + N * 2 ]+E[ i     + N * 2 ];
		FR3 = F[ i     + N * 2 ]+E[ i + 1 + N * 2 ];
		A[ i         ] = A[ i         ] - Z*(FR1-FL1);
		A[ i + N     ] = A[ i + N     ] - Z*(FR2-FL2);
		A[ i + N * 2 ] = A[ i + N * 2 ] - Z*(FR3-FL3);
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if( i == 0 ){
		A[ i         ] =  A[ i + 1         ];
		A[ i + N     ] = -A[ i + 1 + N     ];
		A[ i + N * 2 ] =  A[ i + 1 + N * 2 ];
	}
	if( i == N - 1 ){	
		A[ i         ] =  A[ i - 1         ];
		A[ i + N     ] = -A[ i - 1 + N     ];
		A[ i + N * 2 ] =  A[ i - 1 + N * 2 ];
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if( i < N ){	
		B[ i ] = A[ i ];
		C[ i ] = A[ i + N ]/A[i];
		D[ i ] = ( ( A[ i + N * 2 ]/A[ i ])-0.5*C[ i ]*C[ i ])/CV;
	}	
	barrier( CLK_LOCAL_MEM_FENCE );
}
