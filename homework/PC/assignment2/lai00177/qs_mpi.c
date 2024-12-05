#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <stdint.h>


/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

#define ENABLE_SIZE_CHECK 0
#define SAFTY_FACTOR 4
#define DEBUG 0
#define USE_MEAN 0


int Partition_Single(uint32_t *array, int low, int high, uint32_t pi);
void Quicksort_Single(uint32_t *array, int low, int high);
void Quicksort_Parallel(uint32_t *array, MPI_Comm comm);
void Allocate_Memory(int N);
void Initialize(int N);
void Free_Memory();
void Redistribution(MPI_Comm comm);
void Gather_To_main(MPI_Comm comm);
static void print_numbers(char const * const filename, uint32_t const * const numbers, uint32_t const nnumbers);
static void print_time(double const seconds);
uint32_t Pivot_Parallel(uint32_t p_value, int tag, MPI_Comm comm);
#if ENABLE_SIZE_CHECK
void size_check(int N_check);
#endif


void Swap(uint32_t* array, int i, int j){
    int temp=array[i];
    array[i]=array[j];
    array[j]=temp;
}

static inline double monotonic_seconds()
{
#ifdef __MACH__
  /* OSX */
  static mach_timebase_info_data_t info;
  static double seconds_per_unit;
  if(seconds_per_unit == 0) {
    mach_timebase_info(&info);
    seconds_per_unit = (info.numer / info.denom) / 1e9;
  }
  return seconds_per_unit * mach_absolute_time();
#else
  /* Linux systems */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

int *offset, *n_local, *pivot_list, *n_total;
uint32_t *array, *buffer;
char *output;
int rank, np, np_total, N, NL, n_array, times ,color, old_rank;
double st, et;
#if ENABLE_SIZE_CHECK
int N_MAX;
#endif

int main(int argc, char **argv){
    int i,j, temp;
    double st,et; 
    if(argc<3) {
		printf("Not enough of arguemts\n");
		return 1;
	}
	if(argc>3) {
		printf("Too many arguements\n");
		return 2;
	}
    times=0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    N=atoi(argv[1]);
    output=argv[2];
    temp=(N/np);
    n_array=(rank<(N-temp*np))?temp+1:temp;

#if ENABLE_SIZE_CHECK
    N_MAX=(rank==0)?np*n_array*SAFTY_FACTOR:n_array*SAFTY_FACTOR;
#endif
    //variable here are only for validation
    np_total=np;
    old_rank=rank;
    color=1;

    //if number of points are less than processor, use rank 0 to do quicksort
    if(N<np){
        Allocate_Memory(N);
        if(rank==0){
            st=monotonic_seconds();
            Initialize(N);
            Quicksort_Single(array, 0,N);
            et=monotonic_seconds();
            print_time(et-st);
            print_numbers(output,array,N);

        }
    }else{
        Allocate_Memory(n_array);
        st=monotonic_seconds();
        Initialize(n_array);
        Quicksort_Parallel(array, comm);
#if DEBUG
        for(i=0; i<np; i++){
            if(i==rank){
                printf("at %d, %d, rough final %d(%d): ",times, color, rank, old_rank);
                for(j=0; j<n_array; j++){
                    printf("%d ",array[j]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
#endif    
        Redistribution(comm);
#if DEBUG
        for(i=0; i<np_total; i++){
            if(i==rank){
                printf("final %d: ",rank);
                for(j=0; j<n_array; j++){
                    printf("%d ",array[j]);
                }
                printf("\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
        et=monotonic_seconds();
        if(rank==0){
            print_time(et-st);
        }
        //gather to main and save result
        Gather_To_main(comm);
    }
    Free_Memory();
    MPI_Finalize();
    return 0;
}

void Allocate_Memory(int N){
    size_t size=sizeof(uint32_t);
    if(rank==0){
        array = (uint32_t*)malloc(SAFTY_FACTOR*N*np*size);
        buffer = (uint32_t*)malloc(SAFTY_FACTOR*N*np*size);
    }else{
        array = (uint32_t*)malloc(SAFTY_FACTOR*N*size);
        buffer = (uint32_t*)malloc(SAFTY_FACTOR*N*size);
    }
    size=sizeof(int);
    offset = (int*)malloc(np*size);
    n_local = (int*)malloc(2*(np+1)*size);
    n_total = (int*)malloc(np*size);
}

void Free_Memory(){
    free(array);
    free(n_local);
    free(n_total);
    free(offset);
    free(buffer);
}

void Initialize(int N){
    int i=0;
    srand(time(NULL)+rank);
    for(i=0; i<N; i++){
        array[i]=rand();
    }
}



void Redistribution(MPI_Comm comm){
    int i, j, k;
    int temp, start, red, index, end;
    n_local[0]=0;
    red=0;
    index=0;
    start=0;
    end=start;
#if DEBUG    
    printf("red:%d,  n_array %d\n",rank, n_array);
#endif
    MPI_Scan(&n_array, &temp, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allgather(&temp, 1, MPI_INT, n_local+1 , 1, MPI_INT, comm);
#if DEBUG
    for(i=0; i<np+1; i++){
        printf("red:%d, N_local[%d]=%d\n", rank, i, n_local[i]);
    }
#endif
    temp=(N/np);
    
    for(i=0; i<np; i++){
        
        for(j=0; j<np; j++){
            n_total[j]=0;
        }
        k=(i<(N-temp*np))?temp+1:temp;
        end+=k;
        while(end>n_local[index+1]){
            n_total[index]=n_local[index+1]-start;
            start=n_local[index+1];
            index++;
        }
        n_total[index]=end-start; 
        start=end;
        for(j=0; j<np; j++){
            offset[j]=(j==0)?0:offset[j-1]+n_total[j-1];
        }
        MPI_Gatherv(array+red,n_total[rank],MPI_INT, buffer, n_total, offset, MPI_UINT32_T, i, comm);
        red+=n_total[rank];
    }

    n_array=(rank<(N-temp*np))?temp+1:temp;
    for(i=0; i<n_array; i++){
        array[i]=buffer[i];
    }
}

int Partition_Single(uint32_t *array, int low, int high, uint32_t pivot){
    int j;
    int i = (low-1);
    for(j=low; j<= high; j++){
        if(array[j]<=pivot){
            i++;
            Swap(array, i, j);
        }
    }
    return i;
}

void Quicksort_Parallel(uint32_t *array, MPI_Comm comm){
    int pi;
    int temp, start, red, index, end;
    int npf, npb, nf,nb;
    int tag=1;
    uint32_t ix, pvalue;
    int i,j, k;
    MPI_Comm new_comm;
    if(np==1){
        MPI_Comm_rank(MPI_COMM_WORLD, &temp);
        Quicksort_Single(array, 0, n_array-1);
#if DEBUG
        printf("at %d, %d, %d: original %d with n=%d solved by single processor\n", times, color, rank, temp, n_array);
#endif
    }else{
#if DEBUG
        for(i=0; i<np; i++){
            if(i==rank){
                printf("at %d, %d, original %d(%d): ",times, color, rank, old_rank);
                for(j=0; j<n_array; j++){
                    printf("%d ",array[j]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
        printf("at %d, %d, %d(%d): n_array=%d, np=%d\n",times, color, rank, old_rank, n_array, np);
#endif
        if(n_array){
            pvalue=array[n_array-1];
            tag=1;
        }else{
            pvalue=0;
            tag=0;
        }
        ix=Pivot_Parallel(pvalue, tag, comm);
        pi=Partition_Single(array, 0, n_array-1, ix);
#if DEBUG
        printf("at %d, %d, %d(%d): pi = %d\n",times, color, rank, old_rank, pi);
#endif
        if(n_array){
            nb=pi+1;
            nf=n_array-pi-1;
        }else{
            nb=0;
            nf=0;
        }
#if DEBUG
        for(i=0; i<np; i++){
            if(i==rank){
                printf("at %d, %d, modified %d(%d): ",times, color, rank, old_rank);
                for(j=0; j<n_array; j++){
                    printf("%d ",array[j]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
        printf("at %d, %d, %d(%d): nb=%d, nf=%d\n",times, color, rank, old_rank, nb, nf);
#endif
        //reduction sum of nb, nf for all process
        n_local[0]=0;
        
        MPI_Scan(&nb, &temp, 1, MPI_INT, MPI_SUM, comm);
        MPI_Allgather(&temp, 1, MPI_INT, n_local+1 , 1, MPI_INT, comm);
        n_local[np+1]=0;
        MPI_Scan(&nf, &temp, 1, MPI_INT, MPI_SUM, comm);
        MPI_Allgather(&temp, 1, MPI_INT, n_local+np+2, 1, MPI_INT, comm);
#if DEBUG        
        for(j=0; j<np; j++){
            if(j==rank){
                printf("at %d, %d, %d(%d): nb_local: ",times, color, rank, old_rank);
                for(k=0; k<np+1; k++){
                    printf("%d ",n_local[k]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
        for(j=0; j<np; j++){
            if(j==rank){
                printf("at %d, %d, %d(%d): nf_local: ",times, color, rank, old_rank);
                for(k=0; k<np+1; k++){
                    printf("%d ",n_local[k+np+1]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
        printf("at %d, %d, %d(%d): nb_total=%d, nf_total=%d\n",times, color, rank, old_rank, n_local[np], n_local[2*np+1]);
#endif
        //determine how many process you need
        npb=(int)(np*n_local[np]/(n_local[np]+n_local[2*np+1]));
        npf=np-npb;
        if(npb==0){
            npb++;
            npf--;
        }
        if(npf==0){
            npb--;
            npf++;
        }
#if DEBUG
        printf("at %d, %d, %d(%d): npb=%d, npf=%d, np=%d\n",times, color, rank, old_rank, npb, npf, np);
#endif
        //calculate each process's element
        index=0;
        start=0;
        end=start;
        red=0;
        
#if ENABLE_SIZE_CHECK
        if(rank<npb){
            temp=n_local[np]/npb;
            k=(rank<(n_local[np]-temp*npb))?temp+1:temp;
        }else{
            temp=n_local[2*np+1]/npf;
            k=(rank-npb<(n_local[np]-temp*npf))?temp+1:temp;
        }
        size_check(k);
#endif
        temp=n_local[np]/npb;
        for(i=0; i<npb; i++){
            for(j=0; j<np; j++){
                n_total[j]=0;
            }
            k=(i<(n_local[np]-temp*npb))?temp+1:temp;
            end+=k;
            while(end>n_local[index+1]){
                n_total[index]=n_local[index+1]-start;
                start=n_local[index+1];
                index++;
            }
            n_total[index]=end-start; 
            start=end;
            for(j=0; j<np; j++){
                offset[j]=(j==0)?0:offset[j-1]+n_total[j-1];
            }

            MPI_Gatherv(array+red,n_total[rank],MPI_INT, buffer, n_total, offset, MPI_UINT32_T, i, comm);
            red+=n_total[rank];
#if DEBUG            
            for(j=0; j<np; j++){
                if(j==rank){
                    printf("at %d, %d, %d(%d) to %d nb_total: ", times, color, rank, old_rank, i);
                    for(k=0; k<np; k++){
                        printf("%d ",n_total[k]);
                    }
                    printf("\n");
                }
                MPI_Barrier(comm);
            }
#endif
        }
#if DEBUG
        printf("at %d, %d, %d(%d): backward red=%d\n", times, color, rank, old_rank, red);  
#endif 
        if(rank<npb){
            temp=n_local[np]/npb;
            k=(rank<(n_local[np]-temp*npb))?temp+1:temp;
            n_array=k;
        }
        index=0;
        start=0;
        end=start;
        temp=n_local[2*np+1]/npf;
        for(i=0; i<npf; i++){
            for(j=0; j<np; j++){
                n_total[j]=0;
            }
            k=(i<(n_local[2*np+1]-temp*npf))?temp+1:temp;
            end+=k;
            while(end>n_local[index+np+2]){
                n_total[index]=n_local[index+np+2]-start;
                start=n_local[index+np+2];
                index++;
            }
            n_total[index]=end-start; 
            start=end;

            for(j=0; j<np; j++){
                offset[j]=(j==0)?0:offset[j-1]+n_total[j-1];
            }
            MPI_Gatherv(array+red,n_total[rank],MPI_INT, buffer, n_total, offset, MPI_UINT32_T, i+npb, comm);
            red+=n_total[rank];
#if DEBUG            
            for(j=0; j<np; j++){
                if(j==rank){
                    printf("at %d, %d, %d(%d) to %d nf_total: ", times, color, rank, old_rank, i+npb);
                    for(k=0; k<np; k++){
                        printf("%d ",n_total[k]);
                    }
                    printf("\n");
                }
                MPI_Barrier(comm);
            }
#endif
        }
#if DEBUG
        printf("at %d, %d, %d(%d): forward red=%d\n",times, color,  rank, old_rank, red);
#endif
        if(rank>npb-1){
            temp=n_local[2*np+1]/npf;
            k=((rank-npb)<(n_local[2*np+1]-temp*npf))?temp+1:temp;
            n_array=k;
        }
        for(i=0; i<n_array; i++){
            array[i]=buffer[i];
        }
#if DEBUG        
        for(i=0; i<np; i++){
            if(i==rank){
                printf("at %d, %d, gathered %d(%d): ", times, color, rank, old_rank);
                for(j=0; j<n_array; j++){
                    printf("%d ",array[j]);
                }
                printf("\n");
            }
            MPI_Barrier(comm);
        }
#endif
        //redistribution and create new communitors
        color=(rank<npb)?1:2;
        temp=rank;
        MPI_Comm_split(comm, color, rank, &new_comm);
        MPI_Comm_size(new_comm, &np);
        MPI_Comm_rank(new_comm, &rank);
#if DEBUG        
        printf("at %d, %d, %d(%d): out of %d with n = %d, old =%d\n",times, color, rank, old_rank, np, n_array, temp);
#endif
        //update size and rank
        //np=size of new commucator
        times++;
        Quicksort_Parallel(array, new_comm);
        MPI_Comm_free(&new_comm);
        MPI_Comm_size(comm, &np);
        MPI_Comm_rank(comm, &rank);
    }



}

void Quicksort_Single(uint32_t *array, int low, int high){
    int pi;
    if (low < high){
        pi = Partition_Single(array, low, high, array[high]);
        Quicksort_Single(array, low, pi-1);
        Quicksort_Single(array, pi+1, high);
    }
}

uint32_t Pivot_Parallel(uint32_t p_value, int tag, MPI_Comm comm){
    int i, j, t, temp;
    float sum=0;

    //gather all pivot
    MPI_Allgather(&p_value,1,MPI_INT,buffer,1,MPI_UINT32_T, comm);
    MPI_Allgather(&tag,1,MPI_INT,buffer+np,1,MPI_INT, comm);
#if USE_MEAN

    //find the mean
    for(i=0; i<np; i++){
        sum+=buffer[i];
        temp+=buffer[i+np];
    #if DEBUG       
        if(i==rank){
            printf("at %d, %d, %d(%d), pi: ",times, color, rank, old_rank);
            for(j=0; j<np; j++){
                printf("%d ",buffer[j]);
            }
            printf("\n");
        }
        MPI_Barrier(comm);
    #endif
    }
    sum=sum/temp;
    t=0;
    for(i=1; i<np; i++){
        if(fabs(buffer[t]-sum)>fabs(buffer[i]-sum)){
            t=i;
        }
    }
    
#else
    //find the median
    Quicksort_Single(buffer,0,np-1);
    t=np/2;
#endif
#if DEBUG    
    printf("at %d, %d, %d(%d): sum=%f t=%d, value=%d\n", times, color, rank, old_rank, sum, t, buffer[t]);
#endif
    return buffer[t];
}

void Gather_To_main(MPI_Comm comm){
    int i,j, temp;
    temp=(N/np);
    for(i=0; i<np; i++){
        n_total[i]=(i<(N-temp*np))?temp+1:temp;
        offset[i]=(i==0)?0:offset[i-1]+n_total[i-1];
    }
    MPI_Gatherv(array, n_total[rank], MPI_UINT32_T, buffer, n_total, offset, MPI_UINT32_T, 0, comm);
    if(rank==0){
        print_numbers(output, buffer, N);
    }


}

/**
* @brief Write an array of integers to a file.
*
* @param filename The name of the file to write to.
* @param numbers The array of numbers.
* @param nnumbers How many numbers to write.
*/
static void print_numbers(char const * const filename, uint32_t const * const numbers, uint32_t const nnumbers){
  FILE * fout;
  /* open file */
  if((fout = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write the header */
  fprintf(fout, "%d\n", nnumbers);

  /* write numbers to fout */
  for(uint32_t i = 0; i < nnumbers; ++i) {
    fprintf(fout, "%d\n", numbers[i]);
  }

  fclose(fout);
}

/**
* @brief Output the seconds elapsed while sorting. This excludes input and
*        output time. This should be wallclock time, not CPU time.
*
* @param seconds Seconds spent sorting.
*/
static void print_time(double const seconds){
  printf("Sort Time: %0.04fs\n", seconds);
}

#if ENABLE_SIZE_CHECK
void size_check(int n_check){
    if(N_MAX<n_check){
        uint32_t *temp;
        int i;
        temp=array;
        array=(uint32_t*)malloc(2*N_MAX*sizeof(uint32_t));
        for(i=0; i<n_array; i++){
            array[i]=temp[i];
        }
        free(temp);
        free(buffer);
        buffer=(uint32_t*)malloc(2*N_MAX*sizeof(uint32_t));
        N_MAX=N_MAX*2;
    }
}
#endif