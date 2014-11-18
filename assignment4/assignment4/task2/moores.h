#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <omp.h>

void put_queue(int **vqRef, int vertex, int *qTail, int qSize, int qReaders, omp_lock_t *qLock, int myNum);

int get_queue(int **vqRef, int *qHead, int myNum, int qSize, int *qReaders, omp_lock_t *qLock);

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef, int **dist);

void read_lock(omp_lock_t* lock, int* readers);

void read_unlock(int* readers);

void write_lock(omp_lock_t* lock, int readers);

void write_unlock(omp_lock_t* lock);
