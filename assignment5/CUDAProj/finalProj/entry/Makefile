#gcc -Wall ex1.c -o ex1 -L/home/leon/CUDAProj/finalProj/clamav/lib -lclamav
CLAMPATH = /home/leon/clamav
CUDAPATH = /usr/local/cuda-6.5

NVCC = $(CUDAPATH)/bin/nvcc

CC = gcc
CFLAG = -Wall 
LFLAG = -L$(CLAMPATH)/lib

all: serial

serial: serial.o
	$(CC) -Wall serial.c -o $@ $(LFLAG) -lclamav

clean:
	rm -f *.o serial
