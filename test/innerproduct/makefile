CUDA_INSTALL_PATH :=/usr/local/cuda
CUDA_INCLUDES := -I $(CUDA_INSTALL_PATH)/include
CUDA_LIB := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

all: CPU GPU
	mpiCC main.o cudafunc.o -o test.run $(CUDA_INCLUDES) $(CUDA_LIB)
CPU: 
	mpiCC main.c -c
GPU:
	nvcc cudafunc.cu -c
Clean:
	rm *.dat *.o
