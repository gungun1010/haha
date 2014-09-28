#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

double matmul(int, double*, double*, double*);

int calcSize(int rank, int blockSize);

//initialize the matrix here
void initRowBlk(int sizeA, int sizeB, double** A, double** C);

void initColBlk(int sizeB, double** B);

void matFree(double* A, double* B, double* C);

main(int argc, char **argv) {

    int N, i, run, blockSize, sizeA,sizeB,sizeC;
    double *A, *B, *C, *ArowBlock, *BcolBlock, *CrowBlock;

    int sizes[1];//matrix size
    int p[1];//number of processors 
    double wctime, sparm;
    int rank, procNum, type=99;
    double worktime;
    MPI_Status status;  

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&procNum); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

    /* If I am the master (rank 0) ... */
    if (rank == 0) {
        //sparm = 0; //initialize the workers' work times 
        sizes[0]=4;
        p[0]=2;   
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
        //FIXME
        for (run=0; run<1; run++) {
            
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * N; 
            //init row block A and C based on size
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            
            //send B column block (rank i) to each process with rank i
            for(i=1;i<procNum;i++){

                //variable sizeAB due to different row and col blocks each processor handles
                sizeB = calcSize(i, blockSize); 
                printf("sizeB = %d\n",sizeB); 
                initColBlk(sizeB, &BcolBlock);
                printf("last B = %f\n",BcolBlock[sizeB-1]);
                
                //MPI_Send(dataBuffer, itemNumber, dataType, anotherProcessNum, type, communicator)
                //databuffer: an starting address which holds the data
                //itemNumber: # of items to send starting at the dataBuffer address
                //dataType: the data type of item to send
                //anotherProcessNum: process number of workers
                //type: FIXME,what the????
                //communicator: the pool where the communication happens
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, i, type, MPI_COMM_WORLD);
                
                //recycle B column block memories
                free(BcolBlock);
            }
            //wctime = matmul(blockSize, ArowBlock, BcolBlock, CrowBlock);
/*
            printf ("  %4d     %9.4f\n", N, wctime);
            printf (" serial results: %f\n", C[N*N-1]);
            printf (" ------------------------\n");
*/
        }
    }
    //if im worker for the master
    else{
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
        sizes[0]=200;
        p[0]=2;   
        
        for (run=0; run<1; run++) {
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * N; 
            //init row block A and C based on size
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            
            //variable sizeAB due to different row and col blocks each processor handles
            sizeB = calcSize(rank, blockSize); 
            
            BcolBlock = (double *)calloc(sizeB, sizeof(double));
            
            //FIXME
            MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &status);
            printf("received block B %f at process %d\n", BcolBlock[0],rank);
            //wctime = matmul(blockSize, ArowBlock, BColBlock, CrowBlock);
        }
    }
    
    MPI_Finalize();           

}
