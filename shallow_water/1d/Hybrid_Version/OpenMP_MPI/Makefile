UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
	CCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	CCFLAGS := -fopenmp
endif
ALL: CPU_MPI CPU_CORE 
	mpicc main.o func.o $(CCFLAGS) -lm -o test.run

CPU_MPI:
	mpicc -c $(CCFLAGS) main.c

CPU_CORE:
	mpicc -c $(CCFLAGS) func.c 

clean:
	rm *.o
