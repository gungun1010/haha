#gcc -Wall ex1.c -o ex1 -L/home/leon/CUDAProj/finalProj/clamav/lib -lclamav
CLAMPATH = /home/leon/clamav
CUDAPATH = /usr/local/cuda-6.5

NVCC = $(CUDAPATH)/bin/nvcc

#gcc compiler and clamav library path
CC = gcc
CFLAG = -Wall 
LFLAG = -L$(CLAMPATH)/lib -lclamav

#nvcc compiler and cuda library path
NVCCFLAGS = -I$(CUDAPATH)/include -O3
LFLAGS_NVCC = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

#build for SM50 MAXWELL
GENCODE_SM20 = -g -G -gencode=arch=compute_50,code=\"sm_50,compute_50\"
GENCODE = $(GENCODE_SM20)

.SUFFIXES : .cu .ptx .c

all: serial parallel

serial: serial.o
	$(CC) $(CFLAG) $< -o $@ $(LFLAG)

parallel: parallel.o
	$(NVCC) $(GENCODE) $(LFLAGS_NVCC) -o $@ $< $(LFLAG)


.cu.o:
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

clean:
	rm -f *.o serial
