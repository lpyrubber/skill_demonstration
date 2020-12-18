UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CCFLAGS := -framework OpenCL
else
	CCFLAGS := -lOpenCL
endif

ALL:
	gcc main.c $(CCFLAGS) -o run.exe
