#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(){
    int i,j,k,l,d,N,N2,t;
    int *array, *prefix;
    N=5;
    N2=1;
    while(N2<N){
        N2<<=1;
    }
    N2<<=2;
    array=(int*)malloc(N2*sizeof(int));
    prefix=(int*)malloc((N2+1)*sizeof(int));
    printf("array:\n");
    srand(time(NULL));
    for(i=0; i<N; i++){
        array[i]=rand()%20;
        printf("%d ",array[i]);
    }
    for(i=N; i<N2; i++){
        array[i]=0;
        printf("%d ",array[i]);
    }
    printf("\n");

    //start prefix_sum
    //up sweept
    for(d=2; d<=N2; d<<=1){
        l=d>>1;
        for(k=0; k<N2; k++){
            if((k&(d-1))==0){                
                array[k+d-1]+=array[k+d-l-1];
            }
        }
    }
    printf("after up sweept\n");
    for(i=0; i<N2; i++){
        printf("%d ",array[i]);
    }
    printf("\n");
    //down sweep
    array[N2-1]=0;
    for(d=N2; d>1; d>>=1){
        l=d>>1;
        for(k=0; k<N2; k++){
            if((k&(d-1))==0){
                printf("l=%d,d=%d,k=%d, k-l-1=%d, k-1=%d\n",l,d,k,N2-k-l-1,N2-k-1);
                t=array[N2-1-k-l];
                array[N2-k-1-l]=array[N2-k-1];
                array[N2-k-1]+=t;
            }
        }
        printf("%d: \n",d);
        for(i=0; i<N2; i++){
            printf("%d ",array[i]);
        }
        printf("\n");
    }

    printf("after down sweept\n");
    for(i=0; i<N2; i++){
        printf("%d ",array[i]);
    }
    printf("\n");

    free(array);
    free(prefix);
    return 0;
}