#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void compAndSwap(int *a, int i, int j, int dir);
void bitonicMerge(int *a, int low, int cnt, int dir);
void bitonicSort(int *a,int low, int cnt, int dir);
void sort(int *a, int N, int up);
void swap(int *a, int *b);

int *array;


int main(){
    int i;
    int N, N2;
    N=5;
    N2=1;
    while(N2<N){
        N2*=2;
    }
    printf("N=%d, N2=%d\nbefore\n", N, N2);
    array=(int*)malloc(N2*sizeof(int));
    srand(time(NULL));
    for(i=0; i<N; i++){
        array[i]=rand()%100;
        printf("array[%d]=%d\n",i,array[i]);
    }
    for(i=N; i<N2; i++){
        array[i]=RAND_MAX;
        printf("array[%d]=%d\n",i,array[i]);
    }
    sort(array,N2,1);
    printf("After\n");
    for(i=0; i<N2; i++){
        printf("array[%d]=%d\n",i, array[i]);
    }
    free(array);
    return 0;
}

void compAndSwap(int *a, int i, int j, int dir){
    if (dir==(a[i]>a[j]))
        swap(a+i,a+j);
}
 
/*It recursively sorts a bitonic sequence in ascending order,
  if dir = 1, and in descending order otherwise (means dir=0).
  The sequence to be sorted starts at index position low,
  the parameter cnt is the number of elements to be sorted.*/
void bitonicMerge(int *a, int low, int cnt, int dir){
    if (cnt>1)
    {
        int k = cnt/2;
        for (int i=low; i<low+k; i++)
            compAndSwap(a, i, i+k, dir);
        bitonicMerge(a, low, k, dir);
        bitonicMerge(a, low+k, k, dir);
    }
}

void bitonicSort(int *a,int low, int cnt, int dir){
    if (cnt>1)
    {
        int k = cnt/2;
 
        // sort in ascending order since dir here is 1
        bitonicSort(a, low, k, 1);
 
        // sort in descending order since dir here is 0
        bitonicSort(a, low+k, k, 0);
 
        // Will merge whole sequence in ascending order
        // since dir=1.void sort(int a[], int N, int up)
        bitonicMerge(a,low, cnt, dir);
    }
}
 

void sort(int *a, int N, int up){
    bitonicSort(a,0, N, up);
}

void swap(int *a, int *b){
    int temp;
    temp=*a;
    *a=*b;
    *b=temp;
}