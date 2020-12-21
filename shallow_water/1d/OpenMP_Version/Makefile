UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CCFLAGS := -Xpreprocessor -fopenmp -lomp
else
	CCFLAGS := -fopenmp
endif

all:
	g++ main.c $(CCFLAGS) -o test.run 
