#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include "timing.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
//avoid implicit declaration
void timing(double* wcTime, double* cpuTime);

int calcSize(int rank, int blockSize)
{
    int n;
    int sizeAB;

    n = (rank+1)*blockSize;
    sizeAB = n*(1+n)/2 - (n-blockSize)*(1+n-blockSize)/2;
    return sizeAB;
}

void initRowBlk(int sizeA, int sizeC, double** A, double** C)
{

    int i;
    
    *A = (double *) calloc(sizeA, sizeof(double));//lower triangular mat
    *C = (double *) calloc(sizeC, sizeof(double));//result mat

    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeA; i++) (*A)[i] = 1.0;

}

void initColBlk(int sizeB, double** B)
{

    int i;
    
    *B = (double *) calloc(sizeB, sizeof(double));//upper triangular mat

    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeB; i++) (*B)[i] = 1.0;
    
}

void matFree(double* A, double* B, double* C)
{
    free(A);
    free(B);
    free(C);
}

double matmul(int blockSize, int sizeA, int sizeB, double** A, double** B, double** C) 
{
  int i, j, k;
  int iA, iB, iC;
  int a1,b1;
  double wctime0, wctime1, cputime;

  timing(&wctime0, &cputime);

// This loop computes the matrix-matrix product
  iC = 0;
  a1 = (sizeA + 1 -blockSize)/blockSize-1;
  b1 = (sizeB + 1 -blockSize)/blockSize-1;
  for (i=a1; i<a1+blockSize; i++) {
    //iA = i*(i+1)/2;
    for (j=b1; j<b1+blockSize; j++,iC++) {
      //iB = j*(j+1)/2;
      (*C)[iC] = 0.;
      for (k=0; k<=MIN(i,j); k++) (*C)[iC] += (*A)[k] * (*B)[k]; 
    }
  }
  timing(&wctime1, &cputime);
  return(wctime1 - wctime0);
}
