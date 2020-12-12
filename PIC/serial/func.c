#include "pic.h"

void Compute(){
	int i;
	for( i = 0 ; i < Nt ; ++i ){
		Pushx();
		Hist();
		Efield();
		Pushv();
		Pushx();
	}
}

void Pushx(){
	int i;
	for( i = 0 ; i < N ; ++i ){
		x[ i ] = x[ i ] + v[ i ] * dt * 0.5;
		if( x[ i ] < 0 ){
			x[ i ] = x[ i ] + L;
		}
		if( x[ i ] > L ){
			x[ i ] = x[ i ] - L;
		}
	}
}

void Pushv(){
	int i , j ;
	for( i = 0 ; i < N ; ++i ){
		j = ( int )( x[ i ] / dx );
		v[ i ] = v[ i ] - e[ j ] * dt;
	}
}

void Efield(){
	//CG method
	double b , sum , sum1 , sum2 , sum3 , alfa , beta , temp , error;
	int i , k , up , down;
	char flag = 1;
	for( i = 0 ; i < Nx ; ++i ){
		phi[ i ] = 0;
	}
	for( i = 0 ; i < Nx ; ++i ){
		up = ( i + Nx + 1 ) % Nx;
		down = ( i + Nx - 1 ) % Nx;
		b = ( rho[ i ] - 1 ) * dx * dx;
		r[ i ] = b - phi[ down ] - phi[ up ] + 2*phi[ i ];
		p[ i ] = r [ i ];
	}
	for( i = 0 ; i < Nx ; ++i ){
		up = ( i + Nx + 1 ) % Nx;
                down = ( i + Nx - 1 ) % Nx;
		q[ i ] = r[ up ] + r[ down ] - 2 * r[ i ];
	}
	sum1 = 0;
	sum2 = 0;
	for( i = 0 ; i < Nx ; ++i ){
		sum1 += r[ i ] * r[ i ];
		sum2 += p[ i ] * q[ i ];
	}
	alfa = sum1 / sum2;
	k = 0;

	while( ( k < Itstep ) && ( flag ) ){
		sum = 0;
		for( i = 0 ; i < Nx ; ++i ){
			temp = phi[ i ];
			phi[ i ] += alfa * p[ i ];
			temp -= phi[ i ];
			sum += temp * temp;
			r[ i ] = r[i] - alfa * q[ i ];
		}
		sum3 = 0;
		for( i = 0 ; i < Nx ; ++i ){
			sum3 += r[ i ] * r[ i ];
		}
		beta = sum3 / sum1;
		for( i = 0 ; i < Nx ; ++i ){
			p[ i ] = r[ i ] + beta * p[ i ];
		}
		sum1 = sum3;
		sum2 = 0;
		for( i = 0; i < Nx ; ++i ){
			up = ( i + Nx + 1 ) % Nx;
			down = ( i + Nx -1 ) % Nx;
			q[ i ] = P[ down ] + p[ up ] - 2 * p[ i ];
			sum2 += q[ i ] * p[ i ];
		}
		alfa = sum1 / sum2;
		flag = ( error < Itlimit )? 0 : 1;
		++k;
	}
	for( i = 0 ; i < Nx ; ++i ){
		up = ( i + Nx + 1 ) % Nx;
		down = ( i + Nx - 1 ) % Nx;
		e[ i ] = -0.5 * ( phi[ up ] -phi[down] ) / dx;
	}
}



