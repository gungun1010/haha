#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

double matmul(int, int, int, double**, double**, double**);

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
        
        for (run=0; run<1; run++) {
            
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * blockSize; 
            //init row block A and C based on sizeA and sizeC
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            
            //send B column block (rank i) to each process with rank i
            for(i=1;i<procNum;i++){

                //variable sizeAB due to different row and col blocks each processor handles
                sizeB = calcSize(i, blockSize); 
                initColBlk(sizeB, &BcolBlock);
                
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, i, type, MPI_COMM_WORLD);
                
                //recycle B column block memories
                free(BcolBlock);
            }
            
            //variable sizeAB due to different row and col blocks each processor handles
            //after MPI_Send, work on my own task, the rank=0 
            sizeB = calcSize(rank, blockSize); 
            initColBlk(sizeB, &BcolBlock);
            
            wctime = matmul(blockSize, sizeA, sizeB, &ArowBlock, &BcolBlock, &CrowBlock);
            printf("0: element = %f\n",CrowBlock[0]);
            printf("0: element = %f\n",CrowBlock[1]);
            printf("0: element = %f\n",CrowBlock[2]);
            printf("0: element = %f\n",CrowBlock[3]);
            //MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, type, MPI_COMM_WORLD);
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
        sizes[0]=4;
        p[0]=2;   
        
        for (run=0; run<1; run++) {
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            printf("rank %d, size A = %d\n",rank,sizeA);
            sizeC = blockSize * blockSize; 
            //init row block A and C based on size
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            
            //variable sizeAB due to different row and col blocks each processor handles
            sizeB = calcSize(rank, blockSize); 
            
            BcolBlock = (double *)calloc(sizeB, sizeof(double));
            
            MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE, 0, type, MPI_COMM_WORLD, &status);
            printf("received block B %f at process %d\n", BcolBlock[0],rank);
                       
            wctime = matmul(blockSize,sizeA,sizeB,&ArowBlock, &BcolBlock, &CrowBlock);
            printf("1: element = %f\n",CrowBlock[0]);
            printf("1: element = %f\n",CrowBlock[1]);
            printf("1: element = %f\n",CrowBlock[2]);
            printf("1: element = %f\n",CrowBlock[3]);
            
            /*
            if(rank == procNum-1){
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
            }else{
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, MPI_COMM_WORLD);
            }
            */
            free(BcolBlock);
        }
    }
    
    MPI_Finalize();           

}
