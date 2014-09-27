//#include "timing.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
//avoid implicit declaration
void timing(double* wcTime, double* cpuTime);

double matmul(int N, double* A, double* B, double* C) {

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
