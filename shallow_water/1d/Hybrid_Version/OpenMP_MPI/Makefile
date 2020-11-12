ALL: CPU_MPI CPU_CORE 
	mpicc main.o func.o -o test.run -lm -fopenmp

CPU_MPI:
	mpicc -c -fopenmp main.c

CPU_CORE:
	mpicc -c -fopenmp func.c 

clean:
	rm *.o
