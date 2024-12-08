CUDA_INSTALL_PATH :=/usr/local/cuda-12.6
CUDA_INCLUDES := -I $(CUDA_INSTALL_PATH)/include
CUDA_LIB := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

ALL: GPU CPU
	g++ main.o cudafunc.o -o test.run $(CUDA_INCLUDES) $(CUDA_LIB) -lm
CPU:
	g++ main.c -c
GPU:
	nvcc cudafunc.cu -c
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MPCCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	MPCCFLAGS := -fopenmp
endif

DATA:
	g++ data_generator.c -o datagen.exe
all:SERIAL OPENMP PTHREADS

SERIAL:
	g++ km_serial.cpp -O3 -o km_serial.exe

OPENMP:
	g++ km_openmp.cpp $(MPCCFLAGS) -O3 -o km_openmp.exe
PTHREADS:
	g++ km_pthread.cpp -O3 -o km_pthread.exe
