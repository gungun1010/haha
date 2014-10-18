#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef);

void barrier();

void distributeAB(int rank, int procNum, int blockSize, double** Aref, double** Bref);

void gather(double** rowMatRef, int sizeC, double** Cref);
