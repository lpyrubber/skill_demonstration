#include "main.h"

void Calculate_N2();
void Initial();
void Save_Result();

int *h_a;
int *d_a, *d_b;
int N, N2, N_total;
int BPG; 

int main(){
    int i;
    N=8;
    Calculate_N2();
    Allocate_Memory();
    printf("d_a=%p, h_a=%p\n",d_a, h_a);
    Initial();
    Save_Result();
    Send_To_Device();
    GPU_Compute();
    Send_To_Host();
    Save_Result();
    
    printf("After:\narray:");
    for(i=0; i<N2+1; i++){
        printf("%3d ",h_a[i]);
    }
    printf("\nN_total=%d\n",N_total);

    Free_Memory();
}


void Calculate_N2(){
    N2=0;
    while(N2<N){
        N2+=TPB;
    }
    BPG=((int)((N2+TPB-1)/TPB));
    printf("N=%d, N2=%d\n",N,N2);
}

void Initial(){
    int i;
    N_total=0;
    srand(time(NULL));
    printf("Before\narray:");
    for(i=0; i<N; i++){
        h_a[i]=rand()%20;
        N_total+=h_a[i];
        printf("%3d ",h_a[i]);
    }
    for(i=N; i<N2+1; i++){
        h_a[i]=0;
        printf("%3d ",h_a[i]);
    }
    printf("\n");

}

void Save_Result(){
    int i;
    static int c=0;
    FILE *in;
    if(c==0){
  	    in = fopen("result_1.txt","w");
    }else{
	    in = fopen("result_2.txt","w");
    }
    for (i = 0; i < N2+1; ++i) {
        fprintf(in,"%d\n",  h_a[i]);
    }
    fclose(in);
    c++;
}