#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void printArray(int *arr, int n){
    int i;
    printf("[%d",arr[0]);
    for (i=1; i < n;i++) {
        printf(",%d",arr[i]);
    }
    printf("]\n");
    printf("[%d",arr[n]);
    for (i=1; i < n;i++) {
        printf(",%d",arr[i+n]);
    }
    printf("]\n");
}

void bitonic_sort(int* arr, int n, char flag);

int main() {
    int n, *arr, i,s;
    int N, N2;
    N=8;
    N2=1;
    while(N2<N){
        N2*=2;
    }
    // allocate space and read all the numbers 
    arr = (int *)malloc(2*N2*sizeof(int));
    srand(time(NULL));
    for(i=0; i<N; i++){
        arr[i+N2]=rand()%100;
        arr[i]=i;
        printf("array[%d]=%d\n",i,arr[i]);
    }
    for(i=N; i<N2; i++){
        arr[i+N2]=RAND_MAX;
        arr[i]=i;
        printf("array[%d]=%d\n",i,arr[i]);
    }
    // print array before 
    printArray(arr,N2);

    bitonic_sort(arr,N2,0);
    printArray(arr,N2);
    // do merges
    bitonic_sort(arr,N2,1);

    printArray(arr,N2);
    bitonic_sort(arr,N2,0);
    printArray(arr,N2);
    return 0;
}


void bitonic_sort(int* arr, int n, char flag){
    int k, j;
    printf("mask=%d\n",(n<<1)-1);
    for (int k = 2; k <= n; k <<= 1){
        for (int j = k >> 1; j > 0; j >>= 1){
            for (int l = 0; l < n; l++){
                int i=l;
                int m=i&k;
                if(flag){
                    i+=n;
                }
                int ij = i ^ j;
                if (ij > i){
                    if (m==0){
                        if (arr[i] > arr[ij]){
                            int temp = arr[i];
                            arr[i] = arr[ij];
                            arr[ij] = temp;
                            int i2=(i+n)&((n<<1)-1);
                            int ij2=(ij+n)&((n<<1)-1);
                            temp=arr[i2];
                            arr[i2]=arr[ij2];
                            arr[ij2]=temp;
                        }
                    }else{
                        if (arr[i] < arr[ij]){
                            int temp = arr[i];
                            arr[i] = arr[ij];
                            arr[ij] = temp;
                            int i2=(i+n)&((n<<1)-1);
                            int ij2=(ij+n)&((n<<1)-1);
                            temp=arr[i2];
                            arr[i2]=arr[ij2];
                            arr[ij2]=temp;
                        }
                    }
                }
            }
        }
    }
}

