#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int Partition_Single(int *array, int low, int high, int pi, int *tag);
void Quicksort_Single(int *array, int low, int high);
void Quicksort_parallel(int np, int n_array,int *array);
void Allocate_Memory(int N);
void Initialize(int N);
void Free_Memory();
void Redistribution();
int Pivot_Parallel();
void Swap(int* array, int i, int j){
    int temp=array[i];
    array[i]=array[j];
    array[j]=temp;
}

int *array, *offset, *n_local, *pivot_list, *n_total;
int rank, np, np_total;

int main(int argc, char **argv){
    int i,j; 
    int N=10;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    np_total=np;
    Allocate_Memory(N);
    Initialize(N);
    for(i=0; i<np_total; i++){
        if(i==rank){
            printf("original %d: ",rank);
            for(j=0; j<N; j++){
                printf("%d ",array[j]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    Quicksort_parallel(np, N,array);
    for(i=0; i<np_total; i++){
        if(i==rank){
            printf("modified %d: ",rank);
            for(j=0; j<N; j++){
                printf("%d ",array[j]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    Free_Memory();
    MPI_Finalize();
    return 0;
}

void Allocate_Memory(int N){
    size_t size=sizeof(int);
    array = (int*)malloc(3*N*size);
//    temp = (int*)malloc(N*size);
    offset = (int*)malloc(np*size);
    n_local = (int*)malloc(np*size);
    n_total = (int*)malloc(3*size);
//    printf("Set Memory with %d at %d\n",N, rank);
}

void Free_Memory(){
    free(array);
    free(n_local);
    free(n_total);
    free(offset);
//    free(temp);
//    printf("free Memory at %d\n", rank);
}

void Initialize(int N){
    int i=0;
    srand(time(NULL)+rank);
    for(i=0; i<N; i++){
        array[i]=rand()%10;
    }
}



void Redistribution(){
/*    int i, N;
    printf("%d, ",rank);
    for(i=0; i<rank+1; i++){
        array[i]=rank;
        printf("%d ",array[i]);
    }
    printf("\n");
    for(i=0; i<np; i++){
        n_local[i]=i+1;
        offset[i]=(i==0)?0:offset[i-1]+n_local[i-1];
    }
    MPI_Gatherv(array,rank+1,MPI_INT,temp,n_local,offset,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){
        printf("temp= ");
        for(i=0; i<N; i++){
            printf("%d ", temp[i]);
        }
        printf("\n");
    }
*/
}

int Partition_Single(int *array, int low, int high, int pivot, int *tag){
    int j = low;
    int i = (low-1);
    char flag=1;
    while(j<=high && flag){
        if(array[j]==pivot){
            Swap(array, j, high);
            flag=0;
        }
        j++;
    }
    *tag=0;
    for(j=low; j<= high; j++){
        if(array[j]==pivot){
            *tag=*tag+1;
            
        }else if(array[j]<pivot){
             
            i++;
            Swap(array, i, j);
        }
    }
    i=(i<0)?0:i;
    return i;
}

void Quicksort_parallel(int np, int n_array,int *array){
    int pi;
    int nf, nb;
    int i,j;
    int tag=0;
    if(np==1){
        Quicksort_Single(array, 0, n_array);
    }else{
        i=Pivot_Parallel();
        pi=Partition_Single(array, 0, n_array-1, i, &tag);
        printf("%d, pi = %d, value=%d\n", rank, pi, array[pi]);
        nb=pi+1-tag;
        nf=n_array-1-pi;
        
        printf("%d: nb=%d, np=%d, nf=%d\n",rank, nb, tag, nf);
        //reduction sum of nb, nf for all process
        MPI_Allreduce(&nb, n_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&tag, n_total+1, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nf, n_total+2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        printf("%d: nb_total=%d, np_total=%d, nf_total=%d\n",rank, n_total[0], n_total[1], n_total[2]);
        //determine how many process you nee
        //calculate each process's element
        //redistribution and create new communitors
        //update size and rank
        //np=size of new commucator
        //Quicksort_Parallel(np, new size of aray, array);
    }



}

void Quicksort_Single(int *array, int low, int high){
    int pi,tag;
    if (low < high){
        pi = Partition_Single(array, low, high, array[low], &tag);
        Quicksort_Single(array, low, pi-1);
        Quicksort_Single(array, pi+1, high);
    }
}

int Pivot_Parallel(){
    //find the median of local
    return 5;
}