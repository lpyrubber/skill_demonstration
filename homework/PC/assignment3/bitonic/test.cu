/*
 * Parallel bitonic sort using CUDA.
 * Compile with
 * nvcc -arch=sm_11 bitonic_sort.cu
 * Based on http://www.tools-of-computing.com/tc/CS/Sorts/bitonic_sort.htm
 * License: BSD 3
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/* Every thread gets exactly one value in the unsorted array. */
#define THREADS 512 // 2^9
#define BLOCKS 32768 // 2^15
#define NUM_VALS THREADS*BLOCKS 



void print_elapsed(clock_t start, clock_t stop)
{
  double elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
  printf("Elapsed time: %.3fs\n", elapsed);
}

float random_float()
{
  return (float)rand()/(float)RAND_MAX;
}

void array_print(float *arr, int length) 
{
  int i;
  static int c=0;
  FILE *in;
  if(c==0){
  	in = fopen("result_1.txt","w");
  }else{
	in = fopen("result_2.txt","w");
  }
  for (i = 0; i < length; ++i) {
    fprintf(in,"%1.3f\n",  arr[i]);
  }
  fclose(in);
  c++;
}

void array_fill(float *arr, int length)
{
  srand(time(NULL));
  int i;
  for (i = 0; i < length; ++i) {
    arr[i] = random_float();
  }
}

__global__ void bitonic_sort_step(float *dev_values, int j, int k)
{
  unsigned int i, ixj; /* Sorting partners: i and ixj */
  i = threadIdx.x + blockDim.x * blockIdx.x;
  ixj = i^j;

  /* The threads with the lowest ids sort the array. */
  if ((ixj)>i) {
    if ((i&k)==0) {
      /* Sort ascending */
      if (dev_values[i]>dev_values[ixj]) {
        /* exchange(i,ixj); */
        float temp = dev_values[i];
        dev_values[i] = dev_values[ixj];
        dev_values[ixj] = temp;
      }
    }
    if ((i&k)!=0) {
      /* Sort descending */
      if (dev_values[i]<dev_values[ixj]) {
        /* exchange(i,ixj); */
        float temp = dev_values[i];
        dev_values[i] = dev_values[ixj];
        dev_values[ixj] = temp;
      }
    }
  }
}

/**
 * Inplace bitonic sort using CUDA.
 */
void bitonic_sort(float *values)
{
  float *dev_values;
  size_t size = NUM_VALS * sizeof(float);

  cudaMalloc((void**) &dev_values, size);
  cudaMemcpy(dev_values, values, size, cudaMemcpyHostToDevice);

  dim3 blocks(BLOCKS,1);    /* Number of blocks   */
  dim3 threads(THREADS,1);  /* Number of threads  */

  int j, k;
  /* Major step */
  for (k = 2; k <= NUM_VALS; k <<= 1) {
    /* Minor step */
    for (j=k>>1; j>0; j=j>>1) {
      bitonic_sort_step<<<blocks, threads>>>(dev_values, j, k);
    }
  }
  cudaMemcpy(values, dev_values, size, cudaMemcpyDeviceToHost);
  cudaFree(dev_values);
}


__global__ __forceinline__ void bitonicSort(float* a, float* b);
__global__ __forceinline__ void bitonicBuild(float* a, float* b);
void bitonicBuildRunner(float* a, int size);
void bitonicSortRunner(float* a, int size);

int main(void)
{
  clock_t start, stop;

  float *values = (float*) malloc( NUM_VALS * sizeof(float));
  array_fill(values, NUM_VALS);
  array_print(values, NUM_VALS);
  start = clock();
  bitonicBuildRunner(values, NUM_VALS);
	bitonicSortRunner(values, NUM_VALS);
  stop = clock();
  array_print(values, NUM_VALS);
  print_elapsed(start, stop);
}


void bitonicSortRunner(float* a, int size) {
	// Copy over memory
  float* array;
	int mem = sizeof(float) * size;
	
	int blocks = 1;
	while(blocks != size / 2) {
		// Execution config
        dim3 numBlocks = blocks;
		  dim3 threadsPerBlock = size / blocks / 2;
		
		bitonicSort<<<numBlocks, threadsPerBlock>>>(array, int(size / blocks));
		size *= 2;
	}
}

void bitonicBuildRunner(float* a, int size) {
	int blocks = size / 2;
  int j;
	while(blocks != 1) {
		int i = blocks, blockSize = size * (1 - 1 / blocks);
		while(i != 1) {
			dim3 numBlocks = i, threadsPerBlock = blockSize;
			for(j = 0; j < size; j += blockSize, a++) {
				bitonicBuild<<<numBlocks, threadsPerBlock>>>(a, blockSize, i);
			}
			i /= 2;
		}
		blocks /= 2;
	}
}

/**
 * Applies the bitonic sorting algorithm to each thread. It swaps two
 * elements in the two lists if they are out of place.
 */
__global__ __forceinline__ void bitonicSort(float* a, int blockSize) {
	// First we need to figure out what index each thread will access
    int index = threadIdx.x + blockIdx.x * blockSize;
	atomicMin(&a[i + index], 
		atomicMax(&a[i + blockSize + index], a[i + index]));
	__syncthreads();
}

/**
 * Combines two bitonic sequences together to create a new bitonic sequence.
 * @param a Pointer to the start of the bitonic sequence.
 * @param blockSize The size of each sub-array partition.
 * @param t Determines when to switch between ascending and descending.
 */
__global__ __forceinline__ void bitonicBuild(float* a, int blockSize, int t) {
	int index = threadIdx.x + blockIdx.x * blockSize, x = 0, asc = 1;
	float* b = a + index + (blockSize / 2);
	while(x > index) {
		x += t;
		asc++;
	}
	
	if(asc % 2 == 1) {
		atomicMin(&a[index], atomicMax(&b, a[index]));
	}
	else {
		atomicMax(&a[index], atomicMin(&b, a[index]));
	}
}
		