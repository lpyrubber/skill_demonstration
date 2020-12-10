CUDA_INSTALL_PATH :=/usr/local/cuda-11.0
CUDA_INCLUDES := -I $(CUDA_INSTALL_PATH)/include
CUDA_LIB := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

all: CPU_MPI CPU_CORE 
	mpic++ main.o cudafunc.o -o test.run $(CUDA_INCLUDES) $(CUDA_LIB) -lm

CPU_MPI:
	mpic++ -c main.c

CPU_CORE:
	nvcc -c cudafunc.cu 

clean:
	rm *.o
