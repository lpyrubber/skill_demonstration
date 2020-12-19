ALL: CPU_MPI CPU_CORE 
	mpicc main.o func.o -o test.run -lm

CPU_MPI:
	mpicc -c main.c

CPU_CORE:
	mpicc -c func.c 

clean:
	rm *.o
