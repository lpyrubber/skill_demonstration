#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int Partition_Single(int *array, int low, int high, int pi, int *tag);
void Quicksort_Single(int *array, int low, int high);
void Swap(int* array, int i, int j){
    int temp=array[i];
    array[i]=array[j];
    array[j]=temp;
}

    
int tag=0;
int temp[6]={0,1,0,4,4,6};


int main(int argc, char **argv){
    int i,j; 
    int N=10;
    Quicksort_Single(temp,0,5);
    for(i=0; i<6; i++){
        printf("%d ",temp[i]);
    }
    printf("\n");
    return 0;
}

void Quicksort_Single(int *array, int low, int high){
    int pi, i;
    if (low < high){
        pi = Partition_Single(array, low, high, array[low], &tag);
        printf("n_pivot=%d, pi=%d\n",tag,pi);
        for(i=pi-tag+1; i<pi+1; i++){
            printf("a[%d]=%d\n",i,array[i]);
        }
        printf("nb=%d, np=%d, nf=%d, n=%d\n",pi+1-tag-low,tag,high-pi, high+1-low);
        Quicksort_Single(array, low, pi-tag+1);
        Quicksort_Single(array, pi+1, high);
    }
}

int Partition_Single(int *array, int low, int high, int pivot, int *tag){
    int i = (low-1);
    int k,j;
    char flag=1,temp;
    printf("n=%d, pivot=%d\nbefore: ",high-low+1, pivot);
    for(j=low; j<=high; j++){
        printf("%d ",array[j]);
    }
    printf("\n");
    *tag=0;
    for(j=low; j<= high; j++){
        if(array[j]<=pivot){
            if(array[j]==pivot){
                *tag=*tag+1;
                temp=j;
            }
            i++;
            Swap(array, i, j);
            for(k=0; k<((*tag)-!(temp-j)); k++){
                printf("k=%d, swap %d to %d for %d\n",k, i-k, i-k-1, array[i-k]);
                Swap(array,i-k,i-k-1);
            }
            
        }
    }
    i=(i<0)?0:i;
    printf(" after: ");
    for(j=low; j<=high; j++){
        printf("%d ",array[j]);
    }
    printf("\n");
    return i;
}