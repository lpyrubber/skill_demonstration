CUDA_INSTALL_PATH :=/usr/local/cuda-11.0
CUDA_INCLUDES := -I $(CUDA_INSTALL_PATH)/include
CUDA_LIB := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

ALL: GPU CPU
	g++ main.o cudafunc.o -o test.run $(CUDA_INCLUDES) $(CUDA_LIB) -lm
CPU:
	g++ main.c -c
GPU:
	nvcc cudafunc.cu -c
