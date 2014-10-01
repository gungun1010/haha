#include "mpi.h"
#include "matmul.h"

#define MAT_SIZE 80
#define NUM_PROCESSORS 2

//MPI initialization function
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef){
    
    MPI_Init(argcRef, argvRef);

    MPI_Comm_size(MPI_COMM_WORLD, procNumRef); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, rankRef); // Which process am I?
}

void barrier(){
    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
}

void distributeB(int rank, int procNum, int blockSize, double** Bref){
    int i;
    int sizeB;

    //column block (tag = i) to each process with rank i
    for(i=rank;i<procNum;i++){

        //variable sizeAB due to different row and col blocks each processor handles
        sizeB = calcSize(i, blockSize); 
        initColBlk(sizeB, Bref);
        
        //using target rank as tag
        MPI_Send(*Bref, sizeB, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        
        //recycle B column block memories
        free(*Bref);
    }
}

