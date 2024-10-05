UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	MPCCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	MPCCFLAGS := -fopenmp
endif

DATA:
	g++ data_generator.c -o datagen.exe
SERIAL:
	g++ km_serial.c -O3 -o km_serial

OPENMP:
	g++ km_openmp.c $(MPCCFLAGS) -O2 -o km_openmp
PTHREADS:
	g++ km_pthread.c -O2 -o km_pthread
