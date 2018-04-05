# Makefile for cuda just like the ones for the examples
# Hope it runs on the cuda machine

CUDA_INSTALL_PATH = /opt/cuda
CC = g++
OPTFLAG = -O2 -fomit-frame-pointer -ftree-vectorize -ftree-vectorizer-verbose=0  -funroll-loops
NVCC = ${CUDA_INSTALL_PATH}/bin/nvcc
INCDIR = -I../common/inc/
FLAGS = ${OPTFLAG} -I${CUDA_INSTALL_PATH}/include -Wall -g ${INCDIR}
NVFLAGS = -O2 -I${CUDA_INSTALL_PATH}/include --compiler-options -fno-strict-aliasing --ptxas-options=-v -g ${INCDIR}
BITS = $(shell getconf LONG_BIT)
ifeq (${BITS},64)
	LIBSUFFIX := 64
endif
LFLAGS = -L${CUDA_INSTALL_PATH}/lib${LIBSUFFIX} -lm -lstdc++ -lcudart
CLEAN_FILES = game_of_life cuda_main.o cuda_functions.o


game_of_life: cuda_main.o cuda_functions.o
	${CC} ${LFLAGS} -o $@ $^
	cp $@ ../release

cuda_functions.o: cuda_functions.c
	${CC} -c ${FLAGS} -o $@ $^

cuda_main.o: cuda_main.cu
	${NVCC} ${NVFLAGS} -DUNIX -c $^ -o $@

clean:
	\rm -f $(CLEAN_FILES)
	\rm -f ../release/game_of_life

rebuild: clean game_of_life
