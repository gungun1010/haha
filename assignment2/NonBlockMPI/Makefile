CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

all: nonBlockingMPI
nonBlockingMPI: nonBlockingMPI.o matmul.o ~ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f nonBlockingMPI *.o
