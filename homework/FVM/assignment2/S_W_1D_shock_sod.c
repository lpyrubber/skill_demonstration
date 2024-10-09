#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Gamma 1.4
#define N 1001
#define L 2.0
#define T 0.5
#define dt 0.0001
#define dx ((L)/(N-1)) 
#define N_step (int)(T/dt)
//#define N_step 1

float *U, *Fp, *Fm, *rho, *P, *v;

void Initial();
void Memory_Allocate();
void Compute_Flux();
void Update_U();
void Save_Result();
void Free_Memory();
void Matrix_Multiply(float *a, float *b, float *c, int num);

int main(){
    int i,j,k;

    Memory_Allocate();
    Initial();
    for(i=0; i<N_step; i++){
        Compute_Flux();
        Update_U();
    }
    Save_Result();
    Free_Memory();

}

void Memory_Allocate(){
    size_t num = N*sizeof(float);
    U=(float*)malloc(num*3);
    Fm=(float*)malloc(num*3);
    Fp=(float*)malloc(num*3);
    rho=(float*)malloc(num);
    P=(float*)malloc(num);
    v=(float*)malloc(num);
}

void Initial(){
    int i, j;
    for(i=0; i<N; i++){
        if((-1+dx*i)>0){
            rho[i]=0.125;
            P[i]=0.1;
        }else{
            rho[i]=1.0;
            P[i]=1.0;
        }
        v[i]=0;
        U[i]=rho[i];
        U[i+N]=rho[i]*v[i];
        U[i+2*N]=P[i]/(Gamma-1) + 0.5*rho[i]*v[i]*v[i];
    }
    
}

void Compute_Flux(){
    float S[9],S_inv[9],C[9],C_inv[9],A[9], T_inv[9], LS[9], LT[9],a;
    int i,j,k;
    for(i=0; i<N; i++){
        a=sqrtf(Gamma*P[i]/rho[i]);
        C[0]=1;
        C[1]=0;
        C[2]=0;

        C[3]=-v[i]/rho[i];
        C[4]=1/rho[i];
        C[5]=0;

        C[6]=0.5*(Gamma-1)*v[i]*v[i];
        C[7]=-(Gamma-1)*v[i];
        C[8]=Gamma-1;

        C_inv[0]=1;
        C_inv[1]=0;
        C_inv[2]=0;

        C_inv[3]=v[i];
        C_inv[4]=rho[i];
        C_inv[5]=0;

        C_inv[6]=0.5*v[i]*v[i];
        C_inv[7]=rho[i]*v[i];
        C_inv[8]=1/(Gamma-1);


  

        S_inv[0]=1;
        S_inv[1]=0.5/(a*a);
        S_inv[2]=0.5/(a*a);

        S_inv[3]=0;
        S_inv[4]=0.5/(rho[i]*a);
        S_inv[5]=-0.5/(rho[i]*a);

        S_inv[6]=0;
        S_inv[7]=0.5;
        S_inv[8]=0.5;

        Matrix_Multiply(C_inv, S_inv, T_inv, 3);

        //set up S for Fp
        
        LS[0]=0.5*(v[i]+fabs(v[i]));
        LS[1]=0;
        LS[2]=-0.5*(v[i]+fabs(v[i]))/(a*a);

        LS[3]=0;
        LS[4]=0.5*(v[i]+a+fabs(v[i]+a))*rho[i]*a;
        LS[5]=0.5*(v[i]+a+fabs(v[i]+a));

        LS[6]=0;
        LS[7]=-0.5*(v[i]-a+fabs(v[i]-a))*rho[i]*a;
        LS[8]=0.5*(v[i]-a+fabs(v[i]-a));

        Matrix_Multiply(LS, C, LT, 3);
        Matrix_Multiply(T_inv, LT, A, 3);

        for(j=0;j<3;j++){
            Fp[i+j*N]=0;

            for(k=0;k<3;k++){
                Fp[i+j*N]+=A[k+j*3]*U[i+k*N];
            }
        }
        //set up S for Fm

        LS[0]=0.5*(v[i]-fabs(v[i]));
        LS[1]=0;
        LS[2]=-0.5*(v[i]-fabs(v[i]))/(a*a);

        LS[3]=0;
        LS[4]=0.5*(v[i]+a-fabs(v[i]+a))*rho[i]*a;
        LS[5]=0.5*(v[i]+a-fabs(v[i]+a));

        LS[6]=0;
        LS[7]=-0.5*(v[i]-a-fabs(v[i]-a))*rho[i]*a;
        LS[8]=0.5*(v[i]-a-fabs(v[i]-a));

        Matrix_Multiply(LS, C, LT, 3);
        Matrix_Multiply(T_inv, LT, A, 3);

        for(j=0;j<3;j++){
            Fm[i+j*N]=0;
            for(k=0;k<3;k++){
                Fm[i+j*N]+=A[k+j*3]*U[i+k*N];
            }
        }
    }
}

void Update_U(){
    int i,j,k;
    float FL1, FL2, FL3, FR1, FR2, FR3;
    for(i=1; i<N-1; i++){
        //update U with F

        FL1=Fp[i-1]+Fm[i];
        FR1=Fp[i]+Fm[i+1];
        U[i]=U[i]-(dt/dx)*(FR1-FL1);
        
        FL2=Fp[i-1+N]+Fm[i+N];
        FR2=Fp[i+N]+Fm[i+1+N];
        U[i+N]=U[i+N]-(dt/dx)*(FR2-FL2);

        FL3=Fp[i-1+2*N]+Fm[i+2*N];
        FR3=Fp[i+2*N]+Fm[i+1+2*N];
        U[i+2*N]=U[i+2*N]-(dt/dx)*(FR3-FL3);

        //update rho v P
        rho[i]=U[i];
        v[i]=U[i+N]/U[i];
        P[i]=(Gamma-1)*(U[i+2*N]-0.5*rho[i]*v[i]*v[i]);
    }
}

void Save_Result(){
    FILE *pFile;
    int i,j,k;
    pFile = fopen("data.txt","w");
    for(i=0; i<N; i++){
        fprintf(pFile,"%f %f %f %f\n",-1+dx*i, rho[i], v[i], P[i]);
    }
    fclose(pFile);
}

void Matrix_Multiply(float *a, float *b, float *c, int num){
    int i,j,k;
    for(i=0;i<num;i++){
        for(j=0;j<num;j++){
            c[i+j*num]=0;
            for(k=0;k<num;k++){
                c[i+j*num]+=a[k+j*num]*b[i+k*num];
            }
        }
    }
    
}

void Free_Memory(){
    free(U);
    free(Fm);
    free(Fp);
    free(rho);
    free(P);
    free(v);
}