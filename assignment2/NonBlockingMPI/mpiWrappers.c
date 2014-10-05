#include "mpi.h"
#include "matmul.h"

//MPI initialization function
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef){
    
    MPI_Init(argcRef, argvRef);

    MPI_Comm_size(MPI_COMM_WORLD, procNumRef); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, rankRef); // Which process am I?
}

void barrier(){
    MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
}

void gather(double** rowMatRef, int sizeC, double** Cref){
     barrier(); //wait for everyone to be ready before starting
     MPI_Gather(*rowMatRef,sizeC,MPI_DOUBLE,*Cref,sizeC,MPI_DOUBLE,0,MPI_COMM_WORLD);
     barrier(); //wait for everyone to be ready before starting
}

void distributeAB(int rank, int procNum, int blockSize, double** Aref, double** Bref){
    int i;
    int sizeAB;
    double* bCol;
    double* aRow;
    int bOffset;
    int aOffset;

    bOffset = calcSize(0,blockSize);
    aOffset = bOffset;
    //column block (tag = i) to each process with rank i
    for(i=rank+1;i<procNum;i++){
        //variable sizeAB due to different row and col blocks each processor handles
        sizeAB = calcSize(i, blockSize); 
        initColBlk(sizeAB, &bCol);        
        sliceBlk(&bOffset, sizeAB, &bCol, Bref); 
        
        initRowBlk(sizeAB, &aRow);
        sliceBlk(&aOffset, sizeAB, &aRow, Aref);        
        //using target rank as tag
        MPI_Send(bCol, sizeAB, MPI_DOUBLE, i, i+0xf000, MPI_COMM_WORLD);
        MPI_Send(aRow, sizeAB, MPI_DOUBLE, i, i+0xe000, MPI_COMM_WORLD);
        
        //recycle B column block memories
        free(bCol);
        free(aRow);
    }
}

