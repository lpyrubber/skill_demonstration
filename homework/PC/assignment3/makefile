CUDA_INSTALL_PATH :=/usr/local/cuda-12.6
CUDA_INCLUDES := -I $(CUDA_INSTALL_PATH)/include
CUDA_LIB := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MPCCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	MPCCFLAGS := -fopenmp
endif


ALL: GPU CPU
	g++ km_cuda.o cudafunc.o -o test.run $(CUDA_INCLUDES) $(CUDA_LIB) -lm
CPU:
	g++ km_cuda.cpp -c
GPU:
	nvcc cudafunc.cu -c

DATA:
	g++ data_generator.c -o datagen.exe
all:SERIAL OPENMP PTHREADS

SERIAL:
	g++ km_serial.cpp -O3 -o km_serial.exe

VALID:
	g++ km_openmp.cc $(MPCCFLAGS) -O3 -o validate.exe

OPENMP:
	g++ km_openmp.cpp $(MPCCFLAGS) -O3 -o km_openmp.exe
PTHREADS:
	g++ km_pthread.cpp -O3 -o km_pthread.exe
