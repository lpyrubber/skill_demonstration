UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CCFLAGS := -Xpreprocessor -fopenmp -lomp -mavx
else
	CCFLAGS := -fopenmp -mavx
endif
all:
	g++ main.c $(CCFLAGS) -o test.run
