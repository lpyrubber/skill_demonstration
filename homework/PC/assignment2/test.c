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
int temp[6]={4,1,0,5,4,6};


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
    int pi;
    if (low < high){
        pi = Partition_Single(array, low, high, array[low], &tag);
        printf("n_pivot=%d, pi=%d\n",tag,pi);
//        Quicksort_Single(array, low, pi-1);
//        Quicksort_Single(array, pi+1, high);
    }
}

int Partition_Single(int *array, int low, int high, int pivot, int *tag){
    int j = low;
    int i = (low-1);
    char flag=1;
    while(j<=high && flag){
        if(array[j]==pivot){
            flag=0;
            Swap(array, j, high);
        }
        j++;
    }
    for(j=low; j<= high; j++){
        if(array[j]<=pivot){
            if(array[j]==pivot){
                *tag++;
            }    
            i++;
            Swap(array, i, j);
        }
    }
    return i;
}