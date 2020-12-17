#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <string.h>

#define L       1.0
#define N       200
#define DX      (L/N)
#define DT      (0.01*DX)
#define Z       (DT/DX)
#define NO_STEP 3200
#define R       1.0
#define gamma   ( 7.0 / 5.0 )
#define CV      ( R / ( gamma - 1 ) )
#define CP      ( CV + R )
#define SINGLE  8
#define NS      ((N+SINGLE+1)/SINGLE) 

void Allocate_Memory();
void Initial();
void Compute();
void Save_Result();
void Free_Memory();

float *u1 , *rho;
float *u2 , *v;
float *u3 , *T;
float *fm1 , *fp1;
float *fm2 , *fp2;
float *fm3 , *fp3;
FILE  *pFile;

int main(){
	int i, j;
	pFile = fopen( "data.txt" , "w" );
	Allocate_Memory();
	Initial();
	for(i=0; i<NO_STEP; ++i){
		Compute();
	}
	Save_Result();
	Free_Memory();
	return 0;
}

void Allocate_Memory(){
	int error;
	size_t size , alignment = 32;
	size = ( NS * SINGLE ) * sizeof( float );
	error = posix_memalign( ( void** ) &u1  , alignment , size );
	error = posix_memalign( ( void** ) &rho , alignment , size );
	error = posix_memalign( ( void** ) &u2  , alignment , size );
	error = posix_memalign( ( void** ) &v   , alignment , size );
	error = posix_memalign( ( void** ) &u3  , alignment , size );
	error = posix_memalign( ( void** ) &T   , alignment , size );
	error = posix_memalign( ( void** ) &fm1 , alignment , size );
	error = posix_memalign( ( void** ) &fp1 , alignment , size );
	error = posix_memalign( ( void** ) &fm2 , alignment , size );
	error = posix_memalign( ( void** ) &fp2 , alignment , size );
	error = posix_memalign( ( void** ) &fm3 , alignment , size );
	error = posix_memalign( ( void** ) &fp3 , alignment , size );
}

void Initial(){
	int i;
	for( i = 0 ; i < NS * SINGLE ; ++i ){
		if( i < 0.5 * N  ){
			rho[ i ] = 10.0;	
		}else{
			rho[ i ] = 1.0;
		}
		v[ i ] = 0.0;
		T[ i ] = 1.0;
		u1[ i ] = rho[ i ];
		u2[ i ] = rho[ i ] * v[ i ];
		u3[ i ] = rho[ i ] * ( CV * T[ i ] + 0.5 * v[ i ] * v[ i ] );

	}
}

void Compute(){
	int i;
	__m256 AVX_vel, AVX_a, AVX_F1, AVX_F2, AVX_F3, AVX_Fr;
	__m256 AVX_FL1, AVX_FL2, AVX_FL3, AVX_FR1, AVX_FR2, AVX_FR3; 
	__m256 AVX_t1, AVX_t2, AVX_t3, AVX_t4, AVX_t5, AVX_t6;
	__m256 AVX_cond , AVX_bit;
	//define constant
	__m256 AVX_gamma, AVX_CV , AVX_R , AVX_one, AVX_half, AVX_Z;
	//define array
	__m256 *AVX_u1, *AVX_u2, *AVX_u3, *AVX_rho, *AVX_v, *AVX_T; 
	__m256 *AVX_fp1, *AVX_fp2, *AVX_fp3, *AVX_fm1, *AVX_fm2 , *AVX_fm3;
	AVX_u1  = (__m256*)u1;
	AVX_u2  = (__m256*)u2;
	AVX_u3  = (__m256*)u3;
	AVX_rho = (__m256*)rho;
	AVX_v   = (__m256*)v;
	AVX_T   = (__m256*)T;
	AVX_fp1 = (__m256*)fp1;
	AVX_fp2 = (__m256*)fp2;
	AVX_fp3 = (__m256*)fp3;
	AVX_fm1 = (__m256*)fm1;
	AVX_fm2 = (__m256*)fm2;
	AVX_fm3 = (__m256*)fm3;
	
	AVX_gamma = _mm256_set1_ps( gamma );
	AVX_R     = _mm256_set1_ps( R );
	AVX_CV    = _mm256_set1_ps( CV );
	AVX_one   = _mm256_set1_ps( 1 );
	AVX_half  = _mm256_set1_ps( 0.5 );
	AVX_Z     = _mm256_set1_ps( Z );


	for(i = 0 ; i < NS ; i++){
		// a = sqrt( gamma* R * T[i] );
		AVX_t1 = _mm256_mul_ps( AVX_gamma , AVX_R );
		AVX_t1 = _mm256_mul_ps( AVX_t1 , AVX_T[i] );
		AVX_a  = _mm256_sqrt_ps( AVX_t1 );
		// Fr = v[i]/a;
		AVX_Fr = _mm256_div_ps( AVX_v[i] , AVX_a );
		//if (Fr>1)Fr=1;
		AVX_cond = _mm256_set1_ps( 1 );
		AVX_bit  = _mm256_cmp_ps( AVX_Fr , AVX_cond , _CMP_GT_OS);
		AVX_Fr   = _mm256_blendv_ps( AVX_Fr , AVX_cond , AVX_bit);	
		//if (Fr<-1)Fr=-1;
		AVX_cond = _mm256_set1_ps( -1 );
		AVX_bit  = _mm256_cmp_ps( AVX_Fr , AVX_cond , _CMP_LT_OS);
		AVX_Fr   = _mm256_blendv_ps( AVX_Fr , AVX_cond , AVX_bit);
		//F1 = u1[i]*v[i];
		AVX_F1 = _mm256_mul_ps( AVX_u1[i] , AVX_v[i] );
		//F2 = u1[i]*v[i]*v[i]+rho[i]*R*T[i];
		AVX_t1 = _mm256_mul_ps( AVX_F1 , AVX_v[i] );
		AVX_t2 = _mm256_mul_ps( AVX_rho[i] , AVX_T[i] );
		AVX_t2 = _mm256_mul_ps( AVX_R , AVX_t2 );
		AVX_F2 = _mm256_add_ps( AVX_t1 , AVX_t2);
		//F3 = v[i]*(u3[i]+rho[i]*R*T[i]);
		AVX_t1 = _mm256_add_ps( AVX_u3[i] , AVX_t2 );
		AVX_F3 = _mm256_mul_ps( AVX_v[i] , AVX_t1 ); 
		
		//bit=1-Fr*Fr
		AVX_bit = _mm256_mul_ps( AVX_Fr , AVX_Fr );
		AVX_bit = _mm256_sub_ps( AVX_one , AVX_bit );
		//fp1[i]=0.5*(F1*(Fr+1)+u1[i]*a*bit);
		AVX_t1 = _mm256_add_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F1 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u1[i] , AVX_t2 );
		AVX_fp1[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fp1[i] = _mm256_mul_ps( AVX_half , AVX_fp1[i] );
		//fp2[i]=0.5*(F2*(Fr+1)+u2[i]*a*bit);
		AVX_t1 = _mm256_add_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F2 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u2[i] , AVX_t2 );
		AVX_fp2[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fp2[i] = _mm256_mul_ps( AVX_half , AVX_fp2[i] );
		//fp3[i]=0.5*(F3*(Fr+1)+u3[i]*a*bit);
		AVX_t1 = _mm256_add_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F3 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u3[i] , AVX_t2 );
		AVX_fp3[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fp3[i] = _mm256_mul_ps( AVX_half , AVX_fp3[i] );
		
		//cond=-0.5
		AVX_cond = _mm256_set1_ps( -0.5 );
		//fm1[i]=-0.5*(F1*(Fr-1)+u1[i]*a*bit);
		AVX_t1 = _mm256_sub_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F1 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u1[i] , AVX_t2 );
		AVX_fm1[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fm1[i] = _mm256_mul_ps( AVX_cond , AVX_fm1[i] );
		//fm2[i]=-0.5*(F2*(Fr-1)+u2[i]*a*bit);
		AVX_t1 = _mm256_sub_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F2 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u2[i] , AVX_t2 );
		AVX_fm2[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fm2[i] = _mm256_mul_ps( AVX_cond , AVX_fm2[i] );
		//fm3[i]=-0.5*(F3*(Fr-1)+u3[i]*a*bit);
		AVX_t1 = _mm256_sub_ps( AVX_Fr , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_F3 , AVX_t1 );
		AVX_t2 = _mm256_mul_ps( AVX_a , AVX_bit);
		AVX_t2 = _mm256_mul_ps( AVX_u3[i] , AVX_t2 );
		AVX_fm3[i] = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_fm3[i] = _mm256_mul_ps( AVX_cond , AVX_fm3[i] );
		
	}
	for( i = 0 ; i < NS ; ++i ){
		if( i == 0 ){
			AVX_t1 = _mm256_set_ps( fp1[6], fp1[5], fp1[4], fp1[3]\
					      , fp1[2], fp1[1], fp1[0], 0);
			AVX_t3 = _mm256_set_ps( fp2[6], fp2[5], fp2[4], fp2[3]\
					      , fp2[2], fp2[1], fp2[0], 0);
			AVX_t5 = _mm256_set_ps( fp3[6], fp3[5], fp3[4], fp3[3]\
					      , fp3[2], fp3[1], fp3[0], 0);
		}else{
			AVX_t1 = _mm256_set_ps( fp1[i*8+6], fp1[i*8+5], fp1[i*8+4], fp1[i*8+3]\
					      , fp1[i*8+2], fp1[i*8+1], fp1[i*8+0], fp1[i*8-1]);
			AVX_t3 = _mm256_set_ps( fp2[i*8+6], fp2[i*8+5], fp2[i*8+4], fp2[i*8+3]\
					      , fp2[i*8+2], fp2[i*8+1], fp2[i*8+0], fp2[i*8-1]);
			AVX_t5 = _mm256_set_ps( fp3[i*8+6], fp3[i*8+5], fp3[i*8+4], fp3[i*8+3]\
					      , fp3[i*8+2], fp3[i*8+1], fp3[i*8+0], fp3[i*8-1]);
		}
		if( i == NS - 1){
			AVX_t2 = _mm256_set_ps( 0              , fm1[(NS-1)*8+7], fm1[(NS-1)*8+6], fm1[(NS-1)*8+5]\
					      , fm1[(NS-1)*8+4], fm1[(NS-1)*8+3], fm1[(NS-1)*8+2], fm1[(NS-1)*8+1]);
			AVX_t4 = _mm256_set_ps( 0              , fm2[(NS-1)*8+7], fm2[(NS-1)*8+6], fm2[(NS-1)*8+5]\
					      , fm2[(NS-1)*8+4], fm2[(NS-1)*8+3], fm2[(NS-1)*8+2], fm2[(NS-1)*8+1]);
			AVX_t6 = _mm256_set_ps( 0              , fm3[(NS-1)*8+7], fm3[(NS-1)*8+6], fm3[(NS-1)*8+5]\
					      , fm3[(NS-1)*8+4], fm3[(NS-1)*8+3], fm3[(NS-1)*8+2], fm3[(NS-1)*8+1]);
		}else{

			AVX_t2 = _mm256_set_ps( fm1[i*8+8], fm1[i*8+7], fm1[i*8+6], fm1[i*8+5]\
					      , fm1[i*8+4], fm1[i*8+3], fm1[i*8+2], fm1[i*8+1]);
			AVX_t4 = _mm256_set_ps( fm2[i*8+8], fm2[i*8+7], fm2[i*8+6], fm2[i*8+5]\
					      , fm2[i*8+4], fm2[i*8+3], fm2[i*8+2], fm2[i*8+1]);
			AVX_t6 = _mm256_set_ps( fm3[i*8+8], fm3[i*8+7], fm3[i*8+6], fm3[i*8+5]\
					      , fm3[i*8+4], fm3[i*8+3], fm3[i*8+2], fm3[i*8+1]);
		}
		//FL1=t1+fm1[i]
		AVX_FL1 = _mm256_add_ps( AVX_t1    , AVX_fm1[i] );
		//FR1=fp1[i]+t2
		AVX_FR1 = _mm256_add_ps( AVX_fp1[i] , AVX_t2    );
		//FL2=t3+fm2[i]
		AVX_FL2 = _mm256_add_ps( AVX_t3    , AVX_fm2[i] );
		//FR2=fp2[i]+t4
		AVX_FR2 = _mm256_add_ps( AVX_fp2[i] , AVX_t4    );
		//FL3=t5+fm3[i]
		AVX_FL3 = _mm256_add_ps( AVX_t5    , AVX_fm3[i] );
		//FR3=fp3[i]+t6
		AVX_FR3 = _mm256_add_ps( AVX_fp3[i] , AVX_t6    );
		//u1[i]=u1[i]-Z*(FR1-FL1)
		AVX_t1    = _mm256_sub_ps( AVX_FR1 , AVX_FL1 );
		AVX_t1    = _mm256_mul_ps( AVX_Z , AVX_t1 );
		AVX_u1[i] = _mm256_sub_ps( AVX_u1[i] , AVX_t1);
		//u2[i]=u2[i]-Z*(FR2-FL2)
		AVX_t2    = _mm256_sub_ps( AVX_FR2 , AVX_FL2 );
		AVX_t2    = _mm256_mul_ps( AVX_Z , AVX_t2 );
		AVX_u2[i] = _mm256_sub_ps( AVX_u2[i] , AVX_t2);
		//u3[i]=u3[i]-Z*(FR3-FL3)
		AVX_t3    = _mm256_sub_ps( AVX_FR3 , AVX_FL3 );
		AVX_t3    = _mm256_mul_ps( AVX_Z , AVX_t3 );
		AVX_u3[i] = _mm256_sub_ps( AVX_u3[i] , AVX_t3);
	}
	u1[0]   =  u1[1];
	u2[0]   = -u2[1];
	u3[0]   =  u3[1];
	u1[N+1] =  u1[N];
	u2[N+1] = -u2[N];
	u3[N+1] =  u3[N];
	
	//cond=-0.5
	AVX_cond = _mm256_set1_ps( -0.5 );
	for( i = 0 ; i < NS ; ++i ){
		//rho[i]=u1[i];
		AVX_rho[i] = AVX_u1[i];
		//v[i]=u2[i]/u1[i];
		AVX_v[i] = _mm256_div_ps( AVX_u2[i] , AVX_u1[i] );
		//T[i]=((u3[i]/u1[i])-0.5*v[i]*v[i])/CV;
		AVX_t1   = _mm256_div_ps( AVX_u3[i] , AVX_u1[i] );
		AVX_t2   = _mm256_mul_ps( AVX_v[i] , AVX_v[i] );
		AVX_t2   = _mm256_mul_ps( AVX_cond , AVX_t2 );
		AVX_t2   = _mm256_add_ps( AVX_t1 , AVX_t2 );
		AVX_T[i] = _mm256_div_ps( AVX_t2 , AVX_CV );
	}
}

void Save_Result(){
	int i;
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=1; i<N+1; i++){
			fprintf(pFile, "%e %e %e %e\n",i*DX, rho[i], v[i], T[i]);
		}
	}
}

void Free_Memory(){
	free(fm1);
	free(fp1);
	free(u1);
	free(rho);
	free(fm2);
	free(fp2);
	free(u2);
	free(v);
	free(fm3);
	free(fp3);
	free(u3);
	free(T);
}
