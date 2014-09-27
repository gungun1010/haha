#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

double matmul(int, double*, double*, double*);

//initialize the matrix here
void matInit(int, double* A, double* B, double* C)
{

    int i, sizeAB, sizeC;
    
    sizeAB = N*(N+1)/2; //Only enough space for the nonzero portions of the matrices
    sizeC = N*N; // All of C will be nonzero, in general!

    A = (double *) calloc(sizeAB, sizeof(double));//lower triangular mat
    B = (double *) calloc(sizeAB, sizeof(double));//upper triangular mat
    C = (double *) calloc(sizeC, sizeof(double));//result mat

    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeAB; i++) A[i] = 1.1;
    for (i=0; i<sizeAB; i++) B[i] = 2.1;
}

void matFree(double* A, double* B, double* C)
{
    free(A);
    free(B);
    free(C);
}

int main(int argc, char **argv) {

    /*
    This is the serial program for CPSC424/524 Assignment #2.

    Author: Andrew Sherman, Yale University

    Date: 9/14/2013

    */

    int N, i, run;
    double *A, *B, *C;

    int sizes[4]={100,200,400,800};
    int rowColOffset;
    int p[4]={1,2,4,8};
    double wctime, sparm;
    int rank, size, type=99;
    double worktime;
    MPI_Status status;  

    printf("Matrix multiplication times:\n   N      TIME (secs)\n -----   -------------\n");

    MPI_init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size); // Get no. of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?



    /* If I am the master (rank 0) ... */
    if (rank == 0) {
        sparm = 0; //initialize the workers' work times 
       
        //FIXME
        for (run=0; run<4; run++) {
            N = sizes[run];
            rowColOffset = N/p[run];
             
            matInit(N, A, B, C);

            MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
            
            for(i=0;i<size-1;i++){
                MPI_Send(A[0*rowColOffset], rowColOffset, MPI_DOUBLE, i+1, type, MPI_COMM_WORLD);
                MPI_Send(C[0*rowColOffset], rowColOffset, MPI_DOUBLE, i+1, type, MPI_COMM_WORLD);
            }
            wctime = matmul(N, A, B, C);

            printf ("  %4d     %9.4f\n", N, wctime);
            printf (" serial results: %f\n", C[N*N-1]);
            printf (" ------------------------\n");

        }

        //
    }else{
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
        
        //FIXME
    }
    
    MPI_Finalize();           

}
