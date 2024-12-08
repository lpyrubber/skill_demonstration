#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void merge_up(int *arr, int l, int n) {
    int step=n>>1,i,j,k,temp;
    while (step > 0) {
        for (i=0; i < n; i+=(step<<1)) {
            for (j=i; j < i+step; j++) {
                if (arr[j] > arr[j+step]) {
                    // swap
                    temp = arr[j];
                    arr[j]=arr[j+step];
                    arr[j+step]=temp;
                }
                printf("l=%d, up j=%d, j+step=%d\n",l, j,j+step);
            }
        }
        step >>= 1;
    }
}

void merge_down(int *arr, int l, int n) {
    int step=n>>1,i,j,k,temp;
    while (step > 0) {
        for (i=0; i < n; i+=(step<<1)) {
            for (j=i; j<i+step; j++) {
                if (arr[j] < arr[j+step]) {
                    // swap
                    temp = arr[j];
                    arr[j]=arr[j+step];
                    arr[j+step]=temp;
                }
                printf("l=%d, down j=%d, j+step=%d\n", l,j,j+step);
            }
        }
        step >>= 1;
    }
}

void printArray(int *arr, int n){
    int i;
    printf("[%d",arr[0]);
    for (i=1; i < n;i++) {
        printf(",%d",arr[i]);
    }
    printf("]\n");
}


int main() {
    int n, *arr, i,s;
;   int N, N2;
    N=15;
    N2=1;
    while(N2<N){
        N2*=2;
    }
    // allocate space and read all the numbers 
    arr = (int *)malloc(N2*sizeof(int));
    srand(time(NULL));
    for(i=0; i<N; i++){
        arr[i]=rand()%100;
        printf("array[%d]=%d\n",i,arr[i]);
    }
    for(i=N; i<N2; i++){
        arr[i]=RAND_MAX;
        printf("array[%d]=%d\n",i,arr[i]);
    }
    // print array before 
    printArray(arr,N2);

    // do merges
    for (s=2; s <= N2; s<<=1) {
        for (i=0; i < N2; i+=(s<<1)) {
            merge_up((arr+i),i,s);
            merge_down((arr+i+s),i+s,s);
        }
    }

    printArray(arr,N2);
    return 0;
}