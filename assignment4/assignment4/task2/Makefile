# Make sure to use: module load Langs/Intel/14
CC = icc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -openmp

all: task2
task2:	task2.o moores.o ~ahs3/cpsc424/utils/timing/timing.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o task2 Assignment4RunAll.sh.*
