#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <string.h>
#include <time.h>

#define L       1.0
#define N       20000
#define DX      (L/N)
#define DT      (0.01*DX)
#define Z       (DT/DX)
#define NO_STEP 320000
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
float *f1 , *ML , *MR;
float *f2 , *PL , *PR;
float *f3;
FILE  *pFile;

int main(){
	int i, j;
	clock_t start_t, end_t;
	double total_t;
	start_t=clock();
	pFile = fopen( "data.txt" , "w" );
	Allocate_Memory();
	Initial();
	for(i=0; i<NO_STEP; ++i){
		Compute();
	}
	Save_Result();
	Free_Memory();
	end_t = clock();
	total_t = (double)(end_t-start_t)/ CLOCKS_PER_SEC;
	printf("Total time taken %f sec\n", total_t);
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
	error = posix_memalign( ( void** ) &f1  , alignment , size );
	error = posix_memalign( ( void** ) &f2  , alignment , size );
	error = posix_memalign( ( void** ) &f3  , alignment , size );
	error = posix_memalign( ( void** ) &PL  , alignment , size );
	error = posix_memalign( ( void** ) &PR  , alignment , size );
	error = posix_memalign( ( void** ) &ML  , alignment , size );
	error = posix_memalign( ( void** ) &MR  , alignment , size );
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
	__m256 AVX_a, AVX_M, AVX_Mt, AVX_P;
	__m256 AVX_FL1, AVX_FL2, AVX_FL3, AVX_FR1, AVX_FR2, AVX_FR3; 
	__m256 AVX_t1, AVX_t2, AVX_t3, AVX_t4, AVX_t5, AVX_t6, AVX_t7, AVX_t8, AVX_t9, AVX_t10;
	__m256 AVX_cond , AVX_bit;
	//define constant
	__m256 AVX_gamma, AVX_CV , AVX_R , AVX_one, AVX_half, AVX_Z, AVX_CP, AVX_sign, AVX_quarter;
	//define array
	__m256 *AVX_u1, *AVX_u2, *AVX_u3, *AVX_rho, *AVX_v, *AVX_T; 
	__m256 *AVX_f1, *AVX_f2, *AVX_f3, *AVX_PL, *AVX_PR , *AVX_ML , *AVX_MR;
	AVX_u1  = (__m256*)u1;
	AVX_u2  = (__m256*)u2;
	AVX_u3  = (__m256*)u3;
	AVX_rho = (__m256*)rho;
	AVX_v   = (__m256*)v;
	AVX_T   = (__m256*)T;
	AVX_f1  = (__m256*)f1;
	AVX_f2  = (__m256*)f2;
	AVX_f3  = (__m256*)f3;
	AVX_PL  = (__m256*)PL;
	AVX_PR  = (__m256*)PR;
	AVX_ML  = (__m256*)ML;
	AVX_MR  = (__m256*)MR;
	
	AVX_gamma   = _mm256_set1_ps( gamma );
	AVX_R       = _mm256_set1_ps( R );
	AVX_CV      = _mm256_set1_ps( CV );
	AVX_CP      = _mm256_set1_ps( CP );
	AVX_one     = _mm256_set1_ps( 1 );
	AVX_half    = _mm256_set1_ps( 0.5 );
	AVX_quarter = _mm256_set1_ps( 0.25 );
	AVX_sign    = _mm256_set1_ps( -0.0f );
	AVX_Z       = _mm256_set1_ps( Z );

	for(i = 0 ; i < NS ; i++){
		// a = sqrt( gamma* R * T[i] );
		AVX_t1 = _mm256_mul_ps( AVX_gamma , AVX_R );
		AVX_t1 = _mm256_mul_ps( AVX_t1 , AVX_T[i] );
		AVX_a  = _mm256_sqrt_ps( AVX_t1 );
		// M = v[i]/a;
		AVX_M = _mm256_div_ps( AVX_v[i] , AVX_a );
		// P=rho[i]*R*T[i];
		AVX_P = _mm256_mul_ps( AVX_rho[i] , AVX_R );
		AVX_P = _mm256_mul_ps( AVX_T[i] , AVX_P );
		// Mt=fabs(M);
		AVX_Mt = _mm256_andnot_ps( AVX_sign , AVX_M );
		//( Mt > 1 )?
		AVX_cond = _mm256_set1_ps( 1 );
		AVX_bit  = _mm256_cmp_ps( AVX_Mt , AVX_cond , _CMP_GT_OS);
		
		//t1 = (M+1)*(M+1);
		AVX_t1 = _mm256_add_ps( AVX_M , AVX_one );
		AVX_t1 = _mm256_mul_ps( AVX_t1 , AVX_t1 );
		//t2 = (M-1)*(M-1);
		AVX_t2 = _mm256_sub_ps( AVX_M , AVX_one );
		AVX_t2 = _mm256_mul_ps( AVX_t2 , AVX_t2 );
		//t3 = M+Mt;
		AVX_t3 = _mm256_add_ps( AVX_M , AVX_Mt );
		//t4 = M-Mt;
		AVX_t4 = _mm256_sub_ps( AVX_M , AVX_Mt );
		//cond = 2;
		AVX_cond = _mm256_set1_ps( 2 );
		
		//t5 = 0.5*P*t3/M;
		AVX_t5 = _mm256_div_ps( AVX_t3 , AVX_M );
		AVX_t5 = _mm256_mul_ps( AVX_P , AVX_t5);
		AVX_t5 = _mm256_mul_ps( AVX_half , AVX_t5 );
		//t6 = 0.25*P*t1*(2-M)
		AVX_t6 = _mm256_sub_ps( AVX_cond, AVX_M );
		AVX_t6 = _mm256_mul_ps( AVX_t1, AVX_t6 );
		AVX_t6 = _mm256_mul_ps( AVX_P , AVX_t6 );
		AVX_t6 = _mm256_mul_ps( AVX_quarter , AVX_t6 ); 
		//PL = ( Mt > 1 ) ? t5 : t6
		AVX_PL[i] = _mm256_blendv_ps( AVX_t6 , AVX_t5 , AVX_bit);	
		
		//t5 = 0.5*P*t4/M;
		AVX_t5 = _mm256_div_ps( AVX_t4 , AVX_M );
		AVX_t5 = _mm256_mul_ps( AVX_P , AVX_t5);
		AVX_t5 = _mm256_mul_ps( AVX_half , AVX_t5 );
		//t6 = 0.25*P*t2*(2+M)
		AVX_t6 = _mm256_add_ps( AVX_cond, AVX_M );
		AVX_t6 = _mm256_mul_ps( AVX_t2, AVX_t6 );
		AVX_t6 = _mm256_mul_ps( AVX_P , AVX_t6 );
		AVX_t6 = _mm256_mul_ps( AVX_quarter , AVX_t6 ); 
		//PR = ( Mt > 1 ) ? t5 : t6
		AVX_PR[i] = _mm256_blendv_ps( AVX_t6 , AVX_t5 , AVX_bit);

		//t5 = 0.5*t3;
		AVX_t5 = _mm256_mul_ps( AVX_half  , AVX_t3 );
		//t6 = 0.25*t1;
		AVX_t6 = _mm256_mul_ps( AVX_quarter , AVX_t1 );
		//ML = ( Mt > 1 ) ? t5 : t6
		AVX_ML[i] = _mm256_blendv_ps( AVX_t6 , AVX_t5 , AVX_bit);
		
		//t5 = -0.5*t4;
		AVX_cond = _mm256_set1_ps( -0.5 );
		AVX_t5   = _mm256_mul_ps( AVX_cond  , AVX_t4 );
		//t6 = -0.25*t2;
		AVX_cond = _mm256_set1_ps( -0.25 );
		AVX_t6   = _mm256_mul_ps( AVX_cond , AVX_t2 );
		//ML = ( Mt > 1 ) ? t5 : t6
		AVX_MR[i] = _mm256_blendv_ps( AVX_t6 , AVX_t5 , AVX_bit);
		
		//F1[i]=rho[i]*a;
		AVX_f1[i] = _mm256_mul_ps( AVX_rho[i], AVX_a  );
		//F2[i]=rho[i]*a*v[i];
		AVX_f2[i] = _mm256_mul_ps( AVX_f1[i] , AVX_v[i] );
		//F3[i]=rho[i]*a*(CP*T[i]+0.5*v[i]*v[i]);
		AVX_t1    = _mm256_mul_ps( AVX_v[i] , AVX_v[i] );
		AVX_t1    = _mm256_mul_ps( AVX_half , AVX_t1 );
		AVX_f3[i] = _mm256_mul_ps( AVX_CP , AVX_T[i] );
		AVX_t1    = _mm256_add_ps( AVX_t1 , AVX_f3[i] );
		AVX_f3[i] = _mm256_mul_ps( AVX_f1[i] , AVX_t1 );
		
	}
	for( i = 0 ; i < NS ; ++i ){

		if( i == 0 ){
			AVX_t1 = _mm256_set_ps( f1[6], f1[5], f1[4], f1[3]\
					      , f1[2], f1[1], f1[0], 0);
			AVX_t3 = _mm256_set_ps( f2[6], f2[5], f2[4], f2[3]\
					      , f2[2], f2[1], f2[0], 0);
			AVX_t5 = _mm256_set_ps( f3[6], f3[5], f3[4], f3[3]\
					      , f3[2], f3[1], f3[0], 0);
			AVX_t7 = _mm256_set_ps( ML[6], ML[5], ML[4], ML[3]\
					      , ML[2], ML[1], ML[0], 0);
			AVX_t9 = _mm256_set_ps( PL[6], PL[5], PL[4], PL[3]\
					      , PL[2], PL[1], PL[0], 0);
		}else{
			AVX_t1 = _mm256_set_ps( f1[i*8+6], f1[i*8+5], f1[i*8+4], f1[i*8+3]\
					      , f1[i*8+2], f1[i*8+1], f1[i*8+0], f1[i*8-1]);
			AVX_t3 = _mm256_set_ps( f2[i*8+6], f2[i*8+5], f2[i*8+4], f2[i*8+3]\
					      , f2[i*8+2], f2[i*8+1], f2[i*8+0], f2[i*8-1]);
			AVX_t5 = _mm256_set_ps( f3[i*8+6], f3[i*8+5], f3[i*8+4], f3[i*8+3]\
					      , f3[i*8+2], f3[i*8+1], f3[i*8+0], f3[i*8-1]);
			AVX_t7 = _mm256_set_ps( ML[i*8+6], ML[i*8+5], ML[i*8+4], ML[i*8+3]\
					      , ML[i*8+2], ML[i*8+1], ML[i*8+0], ML[i*8-1]);
			AVX_t9 = _mm256_set_ps( PL[i*8+6], PL[i*8+5], PL[i*8+4], PL[i*8+3]\
					      , PL[i*8+2], PL[i*8+1], PL[i*8+0], PL[i*8-1]);
		}
		if( i == NS - 1){
			AVX_t2  = _mm256_set_ps( 0              , f1[(NS-1)*8+7], f1[(NS-1)*8+6], f1[(NS-1)*8+5]\
					      , f1[(NS-1)*8+4], f1[(NS-1)*8+3], f1[(NS-1)*8+2], f1[(NS-1)*8+1]);
			AVX_t4  = _mm256_set_ps( 0              , f2[(NS-1)*8+7], f2[(NS-1)*8+6], f2[(NS-1)*8+5]\
					      , f2[(NS-1)*8+4], f2[(NS-1)*8+3], f2[(NS-1)*8+2], f2[(NS-1)*8+1]);
			AVX_t6  = _mm256_set_ps( 0              , f3[(NS-1)*8+7], f3[(NS-1)*8+6], f3[(NS-1)*8+5]\
					      , f3[(NS-1)*8+4], f3[(NS-1)*8+3], f3[(NS-1)*8+2], f3[(NS-1)*8+1]);
			AVX_t8  = _mm256_set_ps( 0              , MR[(NS-1)*8+7], MR[(NS-1)*8+6], MR[(NS-1)*8+5]\
					      , MR[(NS-1)*8+4], MR[(NS-1)*8+3], MR[(NS-1)*8+2], MR[(NS-1)*8+1]);
			AVX_t10 = _mm256_set_ps( 0              , PR[(NS-1)*8+7], PR[(NS-1)*8+6], PR[(NS-1)*8+5]\
					      , PR[(NS-1)*8+4], PR[(NS-1)*8+3], PR[(NS-1)*8+2], PR[(NS-1)*8+1]);
		}else{

			AVX_t2  = _mm256_set_ps( f1[i*8+8], f1[i*8+7], f1[i*8+6], f1[i*8+5]\
					      , f1[i*8+4], f1[i*8+3], f1[i*8+2], f1[i*8+1]);
			AVX_t4  = _mm256_set_ps( f2[i*8+8], f2[i*8+7], f2[i*8+6], f2[i*8+5]\
					      , f2[i*8+4], f2[i*8+3], f2[i*8+2], f2[i*8+1]);
			AVX_t6  = _mm256_set_ps( f3[i*8+8], f3[i*8+7], f3[i*8+6], f3[i*8+5]\
					      , f3[i*8+4], f3[i*8+3], f3[i*8+2], f3[i*8+1]);
			AVX_t8  = _mm256_set_ps( MR[i*8+8], MR[i*8+7], MR[i*8+6], MR[i*8+5]\
					      , MR[i*8+4], MR[i*8+3], MR[i*8+2], MR[i*8+1]);
			AVX_t10 = _mm256_set_ps( PR[i*8+8], PR[i*8+7], PR[i*8+6], PR[i*8+5]\
					      , PR[i*8+4], PR[i*8+3], PR[i*8+2], PR[i*8+1]);			  
		}
		//Mt = t7 + MR[i]		
		AVX_Mt = _mm256_add_ps( AVX_t7 , AVX_MR[i] );
		//M = fabs( Mt )
		AVX_M =  _mm256_andnot_ps( AVX_sign , AVX_Mt );
		
		//FL1=Mt*(t1+f1[i])*0.5-M*(f1[i]-t1)*0.5;
		AVX_a   = _mm256_sub_ps( AVX_f1[i] , AVX_t1 );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_FL1 = _mm256_add_ps( AVX_f1[i] , AVX_t1 );
		AVX_FL1 = _mm256_mul_ps( AVX_FL1 , AVX_half );
		AVX_FL1 = _mm256_mul_ps( AVX_FL1 , AVX_Mt );
		AVX_FL1 = _mm256_sub_ps( AVX_FL1 , AVX_a );
		
		//FL2=Mt*(t3+f2[i])*0.5-M*(f2[i]-t3)*0.5+(t9+PR[i]);
		AVX_a   = _mm256_sub_ps( AVX_f2[i] , AVX_t3 );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_P   = _mm256_add_ps( AVX_t9 , AVX_PR[i] );
		AVX_FL2 = _mm256_add_ps( AVX_f2[i] , AVX_t3 );
		AVX_FL2 = _mm256_mul_ps( AVX_FL2 , AVX_half );
		AVX_FL2 = _mm256_mul_ps( AVX_FL2 , AVX_Mt );
		AVX_FL2 = _mm256_sub_ps( AVX_FL2 , AVX_a );
		AVX_FL2 = _mm256_add_ps( AVX_FL2 , AVX_P );
		
		//FL3=Mt*(t5+f3[i])*0.5-M*(f3[i]-t5)*0.5;
		AVX_a   = _mm256_sub_ps( AVX_f3[i] , AVX_t5 );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_FL3 = _mm256_add_ps( AVX_f3[i] , AVX_t5 );
		AVX_FL3 = _mm256_mul_ps( AVX_FL3 , AVX_half );
		AVX_FL3 = _mm256_mul_ps( AVX_FL3 , AVX_Mt );
		AVX_FL3 = _mm256_sub_ps( AVX_FL3 , AVX_a );
		
		//Mt = ML[i] + t8		
		AVX_Mt = _mm256_add_ps( AVX_ML[i] , AVX_t8 );
		//M = fabs( Mt )
		AVX_M =  _mm256_andnot_ps( AVX_sign , AVX_Mt );
		
		//FR1=Mt*(f1[i]+t2)*0.5-fabs(Mt)*(t2-f1[i])*0.5;
		AVX_a   = _mm256_sub_ps( AVX_t2 , AVX_f1[i] );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_FR1 = _mm256_add_ps( AVX_t2 , AVX_f1[i] );
		AVX_FR1 = _mm256_mul_ps( AVX_FR1 , AVX_half );
		AVX_FR1 = _mm256_mul_ps( AVX_FR1 , AVX_Mt );
		AVX_FR1 = _mm256_sub_ps( AVX_FR1 , AVX_a );

		//FR2=Mt*(f2[i]+t4)*0.5-fabs(Mt)*(t4-f2[i])*0.5+(PL[i]+t10);;
		AVX_a   = _mm256_sub_ps( AVX_t4 , AVX_f2[i] );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_P   = _mm256_add_ps( AVX_PL[i] , AVX_t10 );
		AVX_FR2 = _mm256_add_ps( AVX_t4 , AVX_f2[i] );
		AVX_FR2 = _mm256_mul_ps( AVX_FR2 , AVX_half );
		AVX_FR2 = _mm256_mul_ps( AVX_FR2 , AVX_Mt );
		AVX_FR2 = _mm256_sub_ps( AVX_FR2 , AVX_a );
		AVX_FR2 = _mm256_add_ps( AVX_FR2 , AVX_P );

		//FR3=Mt*(f3[i]+t6)*0.5-fabs(Mt)*(t6-f3[i])*0.5;
		AVX_a   = _mm256_sub_ps( AVX_t6 , AVX_f3[i] );
		AVX_a   = _mm256_mul_ps( AVX_half , AVX_a );
		AVX_a   = _mm256_mul_ps( AVX_M , AVX_a );
		AVX_FR3 = _mm256_add_ps( AVX_t6 , AVX_f3[i] );
		AVX_FR3 = _mm256_mul_ps( AVX_FR3 , AVX_half );
		AVX_FR3 = _mm256_mul_ps( AVX_FR3 , AVX_Mt );
		AVX_FR3 = _mm256_sub_ps( AVX_FR3 , AVX_a );
	


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
			fprintf(pFile, "%f %f %f %f\n",i*DX, rho[i], v[i], T[i]);
		}
	}
}

void Free_Memory(){
	free(f1);
	free(f2);
	free(f3);
	free(u1);
	free(u2);
	free(u3);
	free(rho);
	free(v);
	free(T);
	free(ML);
	free(MR);
	free(PL);
	free(PR);
}
