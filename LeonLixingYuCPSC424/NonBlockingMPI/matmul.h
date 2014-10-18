#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

//FIXME, way too many int passing here
double matmul(int rankB, int cCol, int blockSize, int sizeA, int sizeB, double**, double**, double**);

int calcSize(int rank, int blockSize);

void initAB(int sizeAB, double** A, double** B);
//initialize the matrix here
void initRowBlk(int sizeA, double** A);

void initColBlk(int sizeB, double** B);

void sliceBlk(int* offset, int sizeAB, double** AB, double** mat);

void printMat(int size, double** mat);
void debugPrints(double** matRef, int size, int rank);
void updateIndx(int* i, int procNum);

void matFree(double* A, double* B, double* C);

