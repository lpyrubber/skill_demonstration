UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MPCCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	MPCCFLAGS := -fopenmp
endif


ALL: GPU CPU
	nvcc km_cuda.o cudafunc.o -o test.run -lm
CPU:
	g++ km_cuda.cpp -c
GPU:
	nvcc cudafunc.cu -c


VALID:
	g++ km_openmp.cc $(MPCCFLAGS) -O3 -o validate.exe

OPENMP:
	g++ km_openmp.cpp $(MPCCFLAGS) -O3 -o km_openmp.exe
