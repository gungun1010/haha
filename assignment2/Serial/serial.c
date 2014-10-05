#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double matmul(int, double*, double*, double*);

int main(int argc, char **argv) {

    /*
    This is the serial program for CPSC424/524 Assignment #2.

    Author: Andrew Sherman, Yale University

    Date: 9/14/2013

    */

    int N, i, j, k, run;
    double *A, *B, *C;
    int sizeAB, sizeC, iA, iB, iC;

    int sizes[4]={1000,2000,4000,8000};
    
    double wctime;

    printf("Matrix multiplication times:\n   N      TIME (secs)\n -----   -------------\n");
    for (run=0; run<4; run++) {
        N = sizes[run];

        sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrice
        printf("sizeAB = %d\n",sizeAB);
        sizeC = N*N; // All of C will be nonzero, in general!

        A = (double *) calloc(sizeAB, sizeof(double));//lower triangular mat
        B = (double *) calloc(sizeAB, sizeof(double));//upper triangular mat
        C = (double *) calloc(sizeC, sizeof(double));//result mat

        srand(12345); // Use a standard seed value for reproducibility

        // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
         for (i=0; i<sizeAB; i++) A[i] = ((double) rand()/(double)RAND_MAX);
         for (i=0; i<sizeAB; i++) B[i] = ((double) rand()/(double)RAND_MAX);

        wctime = matmul(N, A, B, C);

        printf ("  %4d     %9.4f\n", N, wctime);
        printf ("  last element = %f\n", C[N*N-1]);
        if(N <=16){
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    printf(" %.4f ",C[i*N+j]);
                }
                printf("\n");
            }
        }
        printf (" ------------------------\n");

        free(A);
        free(B);
        free(C);
    }


}
