CUDA_INSTALL_PATH := /usr/local/cuda-10.0
CUDA_LIB          := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

MAIN: CPU_MPI 
	mpicc main.o -o test.run -lm

CPU_MPI:
	mpicc -c main.c

CPU_CORE:
	g++ -c func.c 

clean:
	rm *.o
