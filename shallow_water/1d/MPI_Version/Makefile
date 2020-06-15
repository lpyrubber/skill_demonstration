CUDA_INSTALL_PATH := /usr/local/cuda-10.0
CUDA_LIB          := -L $(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart

ALL: CPU_MPI CPU_CORE 
	mpicc main.o func.o -o test.run -lm

CPU_MPI:
	mpicc -c main.c

CPU_CORE:
	mpicc -c func.c 

clean:
	rm *.o
