#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L     1.0
#define N     200
#define dx    (L/N)
#define dt    0.01*dx
#define Z     (dt/dx)
#define gamma 1.4
#define R     1.0
#define CV    (R/(gamma-1))
#define CP    (CV+R)
#define STEP  3200

float *rho, *u, *T;
float *U, *F;
float *ML, *MR, *PL, *PR;

void Memory_Allocate();
void Initial();
void Compute();
void Save_Result();
void Free_Memory();

int main(){
	int i;
	Memory_Allocate();
	Initial();
	for(i=0; i<STEP; ++i){
		Compute();
	}
	Save_Result();
	Free_Memory();
}

void Memory_Allocate(){
	size_t size=N*sizeof(float);
	rho=(float*)malloc(size);
	u=(float*)malloc(size);
	T=(float*)malloc(size);
	ML=(float*)malloc(size);
	MR=(float*)malloc(size);
	PL=(float*)malloc(size);
	PR=(float*)malloc(size);
	size=3*size;
	U=(float*)malloc(size);
	F=(float*)malloc(size);
}

void Initial(){
	int i;
	for(i=0; i<N; ++i){
		if(i<N*0.5){
			rho[i]=10.0;
			u[i]=0.0;
			T[i]=1.0;
		}else{
			rho[i]=1.0;
			u[i]=0.0;
			T[i]=1.0;
		}
		U[i]=rho[i];
		U[i+N]=rho[i]*u[i];
		U[i+N*2]=rho[i]*(CV*T[i]+0.5*u[i]*u[i]);
	}
}

void Compute(){
	int i;
	float a, M, Mt, P;
	float FL[3], FR[3];
	for(i=0; i<N; ++i){
		a=sqrt(gamma*R*T[i]);
		M=u[i]/a;
		P=rho[i]*R*T[i];
		Mt=fabs(M);
		if(Mt<=1){
			PL[i]= 0.25*P*(M+1)*(M+1)*(2-M);
			PR[i]= 0.25*P*(M-1)*(M-1)*(2+M);
			ML[i]= 0.25*(M+1)*(M+1);
			MR[i]=-0.25*(M-1)*(M-1);
		}else{
			PL[i]= 0.5*P*(M+Mt)/M;
			PR[i]= 0.5*P*(M-Mt)/M;
			ML[i]= 0.5*(M+Mt);
			MR[i]=-0.5*(M-Mt);
		}
		F[i]=rho[i]*a;
		F[i+N]=rho[i]*a*u[i];
		F[i+N*2]=rho[i]*a*(CP*T[i]+0.5*u[i]*u[i]);
//		printf("%f %f %f\n", F[i] , F[i+N] , F[i+N*2] );
	}
	for(i=1; i<N-1; ++i){
		Mt=ML[i-1]+MR[i];
		FL[0]=Mt*(F[i-1    ]+F[i    ])*0.5-fabs(Mt)*(F[i    ]-F[i-1    ])*0.5;
		FL[1]=Mt*(F[i-1+N  ]+F[i+N  ])*0.5-fabs(Mt)*(F[i+N  ]-F[i-1+N  ])*0.5+(PL[i-1]+PR[i]);
		FL[2]=Mt*(F[i-1+N*2]+F[i+N*2])*0.5-fabs(Mt)*(F[i+N*2]-F[i-1+N*2])*0.5;
		Mt=ML[i]+MR[i+1];
		FR[0]=Mt*(F[i    ]+F[i+1    ])*0.5-fabs(Mt)*(F[i+1    ]-F[i    ])*0.5;
		FR[1]=Mt*(F[i+N  ]+F[i+1+N  ])*0.5-fabs(Mt)*(F[i+1+N  ]-F[i+N  ])*0.5+(PL[i]+PR[i+1]);
		FR[2]=Mt*(F[i+N*2]+F[i+1+N*2])*0.5-fabs(Mt)*(F[i+1+N*2]-F[i+N*2])*0.5;
		U[i    ]-=Z*(FR[0]-FL[0]);
		U[i+N  ]-=Z*(FR[1]-FL[1]);
		U[i+N*2]-=Z*(FR[2]-FL[2]);
	}
	U[0  ]= U[1    ];
	U[N  ]=-U[1+N  ];
	U[N*2]= U[1+N*2];
	U[N-1    ]= U[N-2    ];
	U[N-1+N  ]=-U[N-2+N  ];
	U[N-1+N*2]= U[N-2+N*2];
	for(i=0; i<N; ++i){
		rho[i] = U[i];
		u[i] = U[i+N]/U[i];
		T[i] = (U[i+N*2]/U[i]-u[i]*u[i]*0.5)/CV;
	}
}

void Save_Result(){
	FILE *in;
	int i;
	in=fopen("data.txt","w");
	for(i=1; i<N-1; i++){
		fprintf(in, "%f %f %f %f\n", dx*i, rho[i], u[i], T[i]);
	}
	fclose(in);
}

void Free_Memory(){
	free(rho);
	free(u);
	free(T);
	free(ML);
	free(MR);
	free(PL);
	free(PR);
	free(U);
	free(F);
}
