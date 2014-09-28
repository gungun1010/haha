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
    for (i=0; i<sizeA; i++) (*A)[i] = 1.1;

}

void initColBlk(int sizeB, double** B)
{

    int i;
    
    *B = (double *) calloc(sizeB, sizeof(double));//upper triangular mat

    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeB; i++) (*B)[i] = 2.1;
    
}

void matFree(double* A, double* B, double* C)
{
    free(A);
    free(B);
    free(C);
}

double matmul(int N, int rank, int blockSize, double* A, double* B, double* C) {

/*
  This is the serial program for CPSC424/524 Assignment #2.

  Author: Andrew Sherman, Yale University

  Date: 9/14/2014

*/

  int i, j, k;
  int iA, iB, iC;
  double wctime0, wctime1, cputime;

  timing(&wctime0, &cputime);

// This loop computes the matrix-matrix product
  iC = 0;
  if( N!=0 ){
      for (i=0; i<N; i++) {
        iA = i*(i+1)/2;
        for (j=0; j<N; j++,iC++) {
          iB = j*(j+1)/2;
          C[iC] = 0.;
          for (k=0; k<=MIN(i,j); k++) C[iC] += A[iA+k] * B[iB+k]; 
        }
      }
  }else{
  }

  timing(&wctime1, &cputime);
  return(wctime1 - wctime0);
}
