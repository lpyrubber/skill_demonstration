#include <stdlib.h>
#include <stdio.h>

int *a;
int N=10;

int main(){
    int i;
    int *temp;
    a = (int*)malloc(N*sizeof(int));
    for(i=0; i<N; i++){
        a[i]=i;
        printf("%d ",a[i]);
    }
    printf("\n");
    
    temp = a;
    a = (int*)malloc(N*sizeof(int));
    for(i=0; i<N; i++){
        a[i]=temp[i];
    }
    free(temp);
    N=N*2;
    for(i=0; i<N; i++){
        a[i]=i;
        printf("%d ",a[i]);
    }
    printf("\n");
    free(a);
    return 0;
}