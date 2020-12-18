#define gamma 1.4
#define R 1.0
#define CV (R/(gamma-1))
#define CP (CV+R)
__kernel void GPU_Calc(
       		    	int N,
		    	float Z,
		        __global float *A,
		        __global float *B,
		        __global float *C,
		        __global float *D,
		        __global float *E,
		        __global float *F, 
               		__global float *G, 
               		__global float *H, 
                	__global float *I)
{
    int i=get_global_id( 0 );
	float a, M, Mt, P;
	float FL[3], FR[3];
    	if(i<N){
		a=sqrt(gamma*R*D[i]);
		M=C[i]/a;
		P=B[i]*R*D[i];
		Mt=fabs(M);
		if(Mt<=1){
			F[i]= 0.25*P*(M+1)*(M+1)*(2-M);
			G[i]= 0.25*P*(M-1)*(M-1)*(2+M);
			H[i]= 0.25*(M+1)*(M+1);
			I[i]=-0.25*(M-1)*(M-1);
		}else{
			F[i]= 0.5*P*(M+Mt)/M;
			G[i]= 0.5*P*(M-Mt)/M;
			H[i]= 0.5*(M+Mt);
			I[i]=-0.5*(M-Mt);
		}
		E[i]=B[i]*a;
		E[i+N]=B[i]*a*C[i];
		E[i+N*2]=B[i]*a*(CP*D[i]+0.5*C[i]*C[i]);
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if(i>1&&i<N-1){
		Mt=H[i-1]+I[i];
		FL[0]=Mt*(E[i-1    ]+E[i    ])*0.5-fabs(Mt)*(E[i    ]-E[i-1    ])*0.5;
		FL[1]=Mt*(E[i-1+N  ]+E[i+N  ])*0.5-fabs(Mt)*(E[i+N  ]-E[i-1+N  ])*0.5+(F[i-1]+G[i]);
		FL[2]=Mt*(E[i-1+N*2]+E[i+N*2])*0.5-fabs(Mt)*(E[i+N*2]-E[i-1+N*2])*0.5;
		Mt=H[i]+I[i+1];
		FR[0]=Mt*(E[i    ]+E[i+1    ])*0.5-fabs(Mt)*(E[i+1    ]-E[i    ])*0.5;
		FR[1]=Mt*(E[i+N  ]+E[i+1+N  ])*0.5-fabs(Mt)*(E[i+1+N  ]-E[i+N  ])*0.5+(F[i]+G[i+1]);
		FR[2]=Mt*(E[i+N*2]+E[i+1+N*2])*0.5-fabs(Mt)*(E[i+1+N*2]-E[i+N*2])*0.5;
		A[i    ]-=Z*(FR[0]-FL[0]);
		A[i+N  ]-=Z*(FR[1]-FL[1]);
		A[i+N*2]-=Z*(FR[2]-FL[2]);
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if(i==0){
		A[i    ]= A[i+1    ];
		A[i+N  ]=-A[i+1+N  ];
		A[i+N*2]= A[i+1+N*2];
	}
	if(i==N-1){
		A[i    ]= A[i-1    ];
		A[i+N  ]=-A[i-1+N  ];
		A[i+N*2]= A[i-1+N*2];
	}
	barrier( CLK_LOCAL_MEM_FENCE );
	if(i<N){
		B[i] = A[i];
		C[i] = A[i+N]/A[i];
		D[i] = (A[i+N*2]/A[i]-C[i]*C[i]*0.5)/CV;
	}
}                  
