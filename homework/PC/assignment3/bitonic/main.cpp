#include "main.h"

void Calculate_N2();
void Initial();
void Save_Result();

int *h_a;
int *d_a;
int N, N2;
int BPG; 

int main(){
    int i;
    N=462;
    Calculate_N2();
    Allocate_Memory();
    printf("d_a=%p, h_a=%p\n",d_a, h_a);
    Initial();
    Save_Result();
    Send_To_Device();
    GPU_Compute();
    Send_To_Host();
    Save_Result();
/*    
    printf("After:\narray:");
    for(i=0; i<N2; i++){
        printf("%3d ",h_a[i]);
    }
    printf("\n      ");
    for(i=0; i<N2; i++){
        printf("%3d ",h_a[i+N2]);
    }
    printf("\n");
*/
    Free_Memory();
}


void Calculate_N2(){
    N2=1;
    while(N2<N){
        N2=N2*2;
    }
    BPG=((int)((N2+TPB-1)/TPB));
    printf("N=%d, N2=%d\n",N,N2);
}

void Initial(){
    int i;
    srand(time(NULL));
    printf("Before\narray:");
    for(i=0; i<N; i++){
        h_a[i+N2]=rand()%1000;
        h_a[i]=i;
//        printf("%3d ",h_a[i]);
    }
    for(i=N; i<N2; i++){
        h_a[i+N2]=RAND_MAX;
        h_a[i]=i;
//        printf("%3d ",h_a[i]);
    }
/*    printf("\n      ");
    for(i=0; i<N2; i++){
        printf("%3d ",h_a[i+N2]);
    }
    printf("\n");
*/
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
    for (i = 0; i < N2; ++i) {
        fprintf(in,"%d %d\n",  h_a[i], h_a[i+N2]);
    }
    fclose(in);
    c++;
}