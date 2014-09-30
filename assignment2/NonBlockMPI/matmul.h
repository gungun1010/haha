#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#define MAT_SIZE 4
#define NUM_PROCESSORS 4
//FIXME, way too many int passing here
double matmul(int rankB, int cCol, int blockSize, int sizeA, int sizeB, double**, double**, double**);

int calcSize(int rank, int blockSize);

//initialize the matrix here
void initRowBlk(int sizeA, int sizeB, double** A, double** C);

void initColBlk(int sizeB, double** B);

void printMat(int size, double** mat);

void matFree(double* A, double* B, double* C);

