# Make sure to use: module load Langs/Intel/14
CC = icc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

all: task1
task1:	task1.o moores.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o task1
