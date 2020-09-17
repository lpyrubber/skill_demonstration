#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <string.h>

#define L       100.0
#define N       200
#define DX      (L/N)
#define DT      (0.01*DX)
#define Z       (DT/DX)
#define GAP     10
#define NO_STEP 800
#define G       9.81
#define SINGLE  8
#define NS      ((N+SINGLE+1)/SINGLE) 

void Allocate_Memory();
void Initial();
void Compute_Flux();
void Compute_U();
void Update_U();
void Save_Result();
void Free_Memory();

float *u1 , *u1_new;
float *u2 , *u2_new;
float *fm1 , *fp1;
float *fm2 , *fp2;
FILE  *pFile;

int main(){
	int i, j;
	pFile = fopen( "data.txt" , "w" );
	Allocate_Memory();
	Initial();
	for(i=0; i<NO_STEP; ++i){
		Compute_Flux();
//		for(j = 0 ; j < N+2; j++){
//			printf("%e %e %e %e\n",fp1[j], fp2[j], fm1[j], fm2[j]);
//		}
		Compute_U();
		Update_U();
	}
	Save_Result();
	Free_Memory();
	return 0;
}

void Allocate_Memory(){
	int error;
	size_t alignment = 32;

	error = posix_memalign( ( void** ) &u1 ,     alignment , ( NS * SINGLE ) * \
			        sizeof( float ) );
	error = posix_memalign( ( void** ) &u1_new , alignment , ( NS * SINGLE ) * \
		       		sizeof( float ) );
	error = posix_memalign( ( void** ) &u2 ,     alignment , ( NS * SINGLE ) * \
			        sizeof( float ) );
	error = posix_memalign( ( void** ) &u2_new , alignment , ( NS * SINGLE ) * \
		       		sizeof( float ) );
	error = posix_memalign( ( void** ) &fm1 ,    alignment , ( NS * SINGLE ) * \
				sizeof( float ) );
	error = posix_memalign( ( void** ) &fp1 ,    alignment , ( NS * SINGLE ) * \
				sizeof( float ) );	
	error = posix_memalign( ( void** ) &fm2 ,    alignment , ( NS * SINGLE ) * \
				sizeof( float ) );
	error = posix_memalign( ( void** ) &fp2 ,    alignment , ( NS * SINGLE ) * \
				sizeof( float ) );	
}

void Initial(){
	int i;
	for( i = 0 ; i < NS * SINGLE + 2 ; ++i ){
		if( i < 0.5 * N  ){
			u1[ i ] = 10.0;	
		}else{
			u1[ i ] = 1.0;
		}
		u2[ i ] = 0;
		//printf(" %d, %f %f\n", i , u1[i], u2[i]);

	}
}

void Compute_Flux(){
	int i;
	__m256 AVX_vel, AVX_a, AVX_F1, AVX_F2, AVX_Fr;
	__m256 AVX_cond , AVX_bit , AVX_t1, AVX_t2;
	//define constant
	__m256 AVX_G, AVX_one, AVX_half;
	//define array
	__m256 *AVX_u1, *AVX_u2, *AVX_fp1, *AVX_fp2, *AVX_fm1, *AVX_fm2;
	AVX_u1  = (__m256*)u1;
	AVX_u2  = (__m256*)u2;
	AVX_fp1 = (__m256*)fp1;
	AVX_fp2 = (__m256*)fp2;
	AVX_fm1 = (__m256*)fm1;
	AVX_fm2 = (__m256*)fm2;
	AVX_G    = _mm256_set1_ps( G );
	AVX_one  = _mm256_set1_ps( 1 );
	AVX_half = _mm256_set1_ps( 0.5 );
	

	for(i = 0 ; i < NS ; i++){
		//vel = u2[i]/u1[i];
		AVX_vel = _mm256_div_ps( AVX_u2[i] , AVX_u1[i] );
		// a = sqrt(G*u1[i]);
		AVX_t1 = _mm256_mul_ps( AVX_G , AVX_u1[i] );
		AVX_a  = _mm256_sqrt_ps( AVX_t1 );
		// Fr = vel/a;
		AVX_Fr = _mm256_div_ps( AVX_vel , AVX_a );
		//if (Fr>1)Fr=1;
		AVX_cond = _mm256_set1_ps( 1 );
		AVX_bit  = _mm256_cmp_ps( AVX_Fr , AVX_cond , _CMP_GT_OS);
		AVX_Fr   = _mm256_blendv_ps( AVX_Fr , AVX_cond , AVX_bit);	
		//if (Fr<-1)Fr=-1;
		AVX_cond = _mm256_set1_ps( -1 );
		AVX_bit  = _mm256_cmp_ps( AVX_Fr , AVX_cond , _CMP_LT_OS);
		AVX_Fr   = _mm256_blendv_ps( AVX_Fr , AVX_cond , AVX_bit);
		//F1 =u1[i]*vel;
		AVX_F1 = _mm256_mul_ps( AVX_u1[i] , AVX_vel );
		//F2 = F1*vel+0.5*G*u1[i]*u1[i];
		AVX_t1 = _mm256_mul_ps( AVX_F1 , AVX_vel );
		AVX_t2 = _mm256_mul_ps( AVX_u1[i] , AVX_u1[i] );
		AVX_t2 = _mm256_mul_ps( AVX_G , AVX_t2 );
		AVX_t2 = _mm256_mul_ps( AVX_half , AVX_t2 );
		AVX_F2 = _mm256_add_ps( AVX_t1, AVX_t2);
		
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
	}
}

void Compute_U(){
	int i;
	__m256 AVX_FL1, AVX_FL2, AVX_FR1, AVX_FR2, AVX_t1, AVX_t2, AVX_t3, AVX_t4;
	__m256 AVX_Z;
	__m256 *AVX_u1, *AVX_u2, *AVX_u1_new, *AVX_u2_new, *AVX_fp1, *AVX_fp2, *AVX_fm1, *AVX_fm2;
	AVX_u1  = (__m256*)u1;
	AVX_u2  = (__m256*)u2;
	AVX_u1_new = (__m256*)u1_new;
	AVX_u2_new = (__m256*)u2_new;
	AVX_fp1 = (__m256*)fp1;
	AVX_fp2 = (__m256*)fp2;
	AVX_fm1 = (__m256*)fm1;
	AVX_fm2 = (__m256*)fm2;
	AVX_Z = _mm256_set1_ps( Z );
	for( i = 0 ; i < NS ; ++i ){
		if( i == 0 ){
			AVX_t1 = _mm256_set_ps( fp1[6], fp1[5], fp1[4], fp1[3]\
					      , fp1[2], fp1[1], fp1[0], 0);
			AVX_t3 = _mm256_set_ps( fp2[6], fp2[5], fp2[4], fp2[3]\
					      , fp2[2], fp2[1], fp2[0], 0);
		}else{
			AVX_t1 = _mm256_set_ps( fp1[i*8+6], fp1[i*8+5], fp1[i*8+4], fp1[i*8+3]\
					      , fp1[i*8+2], fp1[i*8+1], fp1[i*8+0], fp1[i*8-1]);
			AVX_t3 = _mm256_set_ps( fp2[i*8+6], fp2[i*8+5], fp2[i*8+4], fp2[i*8+3]\
					      , fp2[i*8+2], fp2[i*8+1], fp2[i*8+0], fp2[i*8-1]);
		}
		if( i == NS - 1){
			AVX_t2 = _mm256_set_ps( 0              , fm1[(NS-1)*8+7], fm1[(NS-1)*8+6], fm1[(NS-1)*8+5]\
					      , fm1[(NS-1)*8+4], fm1[(NS-1)*8+3], fm1[(NS-1)*8+2], fm1[(NS-1)*8+1]);
			AVX_t4 = _mm256_set_ps( 0              , fm2[(NS-1)*8+7], fm2[(NS-1)*8+6], fm2[(NS-1)*8+5]\
					      , fm2[(NS-1)*8+4], fm2[(NS-1)*8+3], fm2[(NS-1)*8+2], fm2[(NS-1)*8+1]);
		}else{

			AVX_t2 = _mm256_set_ps( fm1[i*8+8], fm1[i*8+7], fm1[i*8+6], fm1[i*8+5]\
					      , fm1[i*8+4], fm1[i*8+3], fm1[i*8+2], fm1[i*8+1]);
			AVX_t4 = _mm256_set_ps( fm2[i*8+8], fm2[i*8+7], fm2[i*8+6], fm2[i*8+5]\
					      , fm2[i*8+4], fm2[i*8+3], fm2[i*8+2], fm2[i*8+1]);
		}
		//FL1=t1+fm1[i]
		AVX_FL1 = _mm256_add_ps( AVX_t1    , AVX_fm1[i] );
		//FR1=fp1[i]+t2
		AVX_FR1 = _mm256_add_ps( AVX_fp1[i] , AVX_t2    );
		//FL2=t3+fm2[i]
		AVX_FL2 = _mm256_add_ps( AVX_t3    , AVX_fm2[i] );
		//FR2=fp2[i]+t4
		AVX_FR2 = _mm256_add_ps( AVX_fp2[i] , AVX_t4    );
		//u1_new[i]=u1[i]-Z*(FR1-FL1)
		AVX_t1        = _mm256_sub_ps( AVX_FR1 , AVX_FL1 );
		AVX_t1        = _mm256_mul_ps( AVX_Z , AVX_t1 );
		AVX_u1_new[i] = _mm256_sub_ps( AVX_u1[i] , AVX_t1);
		//u2_new[i]=u2[i]-Z*(FR2-FL2)
		AVX_t2        = _mm256_sub_ps( AVX_FR2 , AVX_FL2 );
		AVX_t2        = _mm256_mul_ps( AVX_Z , AVX_t2 );
		AVX_u2_new[i] = _mm256_sub_ps( AVX_u2[i] , AVX_t2);
	}
	u1_new[0]   =  u1_new[1];
	u2_new[0]   = -u2_new[1];
	u1_new[N+1] =  u1_new[N];
	u2_new[N+1] = -u2_new[N];
	
}

void Update_U(){
	int i;
	for( i = 0 ; i < N + 2 ; ++i ){
		u1[i] = u1_new[i];
		u2[i] = u2_new[i];
	}
}

void Save_Result(){
	int i;
	if (pFile == NULL){
		printf("File open failed\n");
	}else{
		for(i=1; i<N+1; i++){
			fprintf(pFile, "%e %e %e\n",i*DX, u1[i], u2[i]);
//			printf( "%e %e %e\n", i*DX , u1[i] , u2[i]);
		}
		fprintf(pFile,"\n");
	}
}

void Free_Memory(){
	free(fm1);
	free(fp1);
	free(u1);
	free(u1_new);
	free(fm2);
	free(fp2);
	free(u2);
	free(u2_new);
}
