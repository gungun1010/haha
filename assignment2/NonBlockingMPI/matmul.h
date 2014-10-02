#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

//FIXME, way too many int passing here
double matmul(int rankB, int cCol, int blockSize, int sizeA, int sizeB, double**, double**, double**);

int calcSize(int rank, int blockSize);

//initialize the matrix here
void initRowBlk(int sizeA, int sizeB, double** A, double** C);

void initColBlk(int sizeB, double** B);


void initAnC(int rank, int blockSize, int N, int* sizeAref,int* sizeCref, double** A, double** C);

void printMat(int size, double** mat);
void debugPrints(double** matRef, int size, int rank);
void updateIndx(int* i, int procNum);

void matFree(double* A, double* B, double* C);

