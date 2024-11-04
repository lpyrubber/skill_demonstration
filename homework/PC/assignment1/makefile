
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
