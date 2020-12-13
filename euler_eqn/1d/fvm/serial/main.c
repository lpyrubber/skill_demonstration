#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N       200
#define L       1.0
#define DX      ( N / L )
#define DT      ( 0.01*DX )
#define Z       ( DT / DX )
#define NO_STEP 1
#define R       1.0
#define gamma   ( 7.0 / 5.0 )
#define CV      ( R / ( gamma - 1 ) )
#define CP      ( CV + R )

void Memory_Allocate();
void Initial();
void Compute();
void Save_Result();
void Free_Memory();

float *u , *fm , *fp;
float *rho , *v , *T;

int main(){
	int i;
	Memory_Allocate();
	Initial();
	for(i=0; i<NO_STEP; i++){
		Compute();
	}
	printf("%f %f %f %f\n", fm[0], fm[1], fp[1], fp[2]);
	Save_Result();
	Free_Memory();
}

void Memory_Allocate(){
	size_t size;
	size = ( N + 2 ) * sizeof( float );
	u   = ( float* )malloc( 3 * size );
	fm  = ( float* )malloc( 3 * size );
	fp  = ( float* )malloc( 3 * size );
	rho = ( float* )malloc( size );
	v   = ( float* )malloc( size );
	T   = ( float* )malloc( size );
}

void Initial(){
	int i;
	for( i = 0 ; i < N + 2 ; ++i ){
		if( i < 0.5 * N ){
			rho[ i ] = 10.0;
			v[ i ]   = 0.0;
			T[ i ]   = 1.0;
		}else{
			rho[ i ] = 1.0;
			v[ i ]   = 0.0;
			T[ i ]   = 1.0;
		}
		u[ i                 ] = rho[ i ];
		u[ i + ( N + 2 )     ] = rho[ i ] * v[ i ];
		u[ i + ( N + 2 ) * 2 ] = rho[ i ] * ( CV * T[ i ] \
				       + 0.5 * v[ i ] * v[ i ] ); 
	}
}

void Compute(){
	int i;
	float vel, a;
	float F1, F2, F3;
	float FL1, FL2, FL3, FR1, FR2, FR3;
	float Fr;
	for( i = 0; i < N + 2 ; ++i ){
		a = sqrt( gamma * R * T[ i ] );
		Fr = v[ i ] / a;
		if ( Fr > 1  ) Fr =  1;
		if ( Fr < -1 ) Fr = -1;
		F1 = u[ i ] * v[ i ];
		F2 = u[ i ] * v[ i ] * v[ i ] +  rho[ i ] * R * T[ i ];
		F3 = v[ i ] * ( u[ i + ( N + 2 ) * 2 ] + rho[ i ] * R * T[ i ] );

		fp[i] = 0.5*(F1*(Fr+1)+u[i]*a*(1-Fr*Fr));
		fp[i+N+2] = 0.5*(F2*(Fr+1)+u[i+N+2]*a*(1-Fr*Fr));
		fp[i+(N+2)*2] = 0.5*(F3*(Fr+1)+u[i+(N+2)*2]*a*(1-Fr*Fr));
		fm[i] = -0.5*(F1*(Fr-1)+u[i]*a*(1-Fr*Fr));
		fm[i+N+2] = -0.5*(F2*(Fr-1)+u[i+N+2]*a*(1-Fr*Fr));
		fm[i+(N+2)*2] = -0.5*(F3*(Fr-1)+u[i+(N+2)*2]*a*(1-Fr*Fr));
	}
	for(i=1; i<N+1; i++){
		FL1 = fp[i-1]+fm[i];
		FR1 = fp[i]+fm[i+1];
		FL2 = fp[i-1+N+2]+fm[i+N+2];
		FR2 = fp[i+N+2]+fm[i+1+N+2];
		FL3 = fp[i-1+(N+2)*2]+fm[i+(N+2)*2];
		FR3 = fp[i+(N+2)*2]+fm[i+1+(N+2)*2];
		u[i] = u[i] - Z*(FR1-FL1);
		u[i+N+2] = u[i+N+2] - Z*(FR2-FL2);
		u[i+(N+2)*2] = u[i+(N+2)*2] - Z*(FR3-FL3);
	}
	
	u[ 0         ] =  u[1];
	u[ 0+N+2     ] = -u[1+N+2];
	u[ 0+(N+2)*2 ] =  u[1+(N+2)*2];
		
	u[ N+1         ] =  u[N];
	u[ N+1+N+2     ] = -u[N+N+2];
	u[ N+1+(N+2)*2 ] =  u[N+(N+2)*2];

	for(i=0; i<N+2; i++){
		rho[i]=u[i];
		v[i]=u[i+N+2]/u[i];
		T[i]=((u[i+(N+2)*2]/u[i])-0.5*v[i]*v[i])/CV;
	}
	
}
void Save_Result(){
	FILE *pFile;
	int i;
	pFile=fopen("data.txt","w");
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=0; i<N+2; i++){
			fprintf(pFile, "%e %e %e %e\n",i*DX, u[i], v[i], T[i]);
		}
		fprintf(pFile,"\n");
	}
	fclose(pFile);
}

void Free_Memory(){
	free(fm);
	free(fp);
	free(u);
	free(rho);
	free(v);
	free(T);
}
