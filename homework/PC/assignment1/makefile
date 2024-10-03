UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MPCCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	MPCCFLAGS := -fopenmp
endif

DATA:
	g++ data_generator.c -o datagen.exe

OPENMP:
	g++ km_openmp.c $(MPCCFLAGS) -o test.run
