#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
#include <omp.h>

void put_queue(int **vqRef, int vertex, int *qSize, omp_lock_t *lock, int *readers);

int get_queue(int **vqRef, int *qSize, omp_lock_t *lock, int *readers);

void read_lock(omp_lock_t* lock, int* readers);

void read_unlock(int* readers);

void write_lock(omp_lock_t* lock, int readers);

void write_unlock(omp_lock_t* lock);

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef, int **dist);
