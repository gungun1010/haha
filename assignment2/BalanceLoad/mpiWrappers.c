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

void distributeB(int rank, int procNum, int blockSize, double** Bref, MPI_Request* requestRef){
    int i;
    int sizeB;
    
    //column block (tag = i) to each process with rank i
    for(i=rank+1;i<procNum;i++){

        //variable sizeAB due to different row and col blocks each processor handles
        sizeB = calcSize(i, blockSize); 
        initColBlk(sizeB, Bref);
        
        //using target rank as tag
        MPI_Send(*Bref, sizeB, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        
        //recycle B column block memories
        free(*Bref);
    }
}

