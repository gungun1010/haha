#include "mpi.h"
#include "matmul.h"
#include "mpiWrappers.h"
//this include is for easy change of variable in the testEntry folder
#include "../testEntry/variables.h"

main(int argc, char **argv) {
    int N, i, run, blockSize, sizeA,sizeB,sizeC,sizeT;
    double *A,*B,*C,*ArowBlock,*BtempBlock,*BcolBlock,*CrowBlock;
    int sizes[1];//matrix size
    int p[1];//number of processors 
    double wctime, sparm;
    int rank, procNum;
    int srcRank, destRank,sendTag;
    double worktime;
    MPI_Status status;  
    MPI_Request request;
    MPI_Request requestBackup;
    init(&argc, &argv, &procNum, &rank);

    /* If I am the master (rank 0) ... */
    if (rank == 0) {
        //sparm = 0; //initialize the workers' work times 
        barrier(); //wait for everyone to be ready before starting timer
        wctime = 0.0;
        sizes[0]=MAT_SIZE;
        p[0]=NUM_PROCESSORS;
        
        for (run=0; run<1; run++) {
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows/cols per block A/:q
            
            //initialize A,B, C block
            sizeT = calcSize(rank,N);
            initAB(sizeT, &A, &B);
            
            //distribute B Col block to each of the workers
            distributeAB(rank, procNum, blockSize, &A, &B);
            
            sizeA = calcSize(rank, blockSize);
            initRowBlk(sizeA, &ArowBlock);
            memcpy(ArowBlock, A, sizeA*sizeof(double)); 
            sizeB = sizeA;
            initColBlk(sizeB, &BcolBlock);
            memcpy(BcolBlock, B, sizeB*sizeof(double)); 
            sizeC=blockSize*N;
            initRowBlk(sizeC, &CrowBlock);
            //after MPI_Send, work on my own task, rank=0
            
            //printf("%d send msg",rank);  
            if(procNum > 1){ 
                MPI_Isend(BcolBlock, sizeB, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD, &request);
            }
            //after MPI_Isend, work on my own task, rank=0
            //parallized matmul and send for large buffer size
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
            
            if(procNum > 1){ 
                MPI_Wait(&request, &status); // wait for send and matmul to finish
            }
            free(BcolBlock); //up to here, B column block buffer has done its job
            
            for(srcRank=procNum-1; srcRank>0; srcRank--){                 
                sizeT = calcSize(srcRank, blockSize); 
                BtempBlock = (double *)calloc(sizeT,sizeof(double));
                
                //recall: I use rank # as my tag, but master only receives message from last worker
                MPI_Irecv(BtempBlock, sizeT, MPI_DOUBLE,procNum-1,srcRank,MPI_COMM_WORLD,&request);
                //must wait for BtempBlock buffer to be filled before matmul can use it
                MPI_Wait(&request, &status);

                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                //we dont send this B back to the original owner
                //the Isend() is pipelined with matmul so with larger load, it is parallel
                if(srcRank != rank+1){
                    MPI_Isend(BtempBlock, sizeT, MPI_DOUBLE, rank+1, srcRank, MPI_COMM_WORLD, &request);
                }//NOTE the differences here, im passing Temp B buffer instead of B Column buff
                
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                MPI_Wait(&request, &status);
                free(BtempBlock);
            }
            
            //collective MPI 
            C = (double *)malloc(N*N*sizeof(double));
            gather(&CrowBlock, sizeC, &C); 
            printf("proc %d: N = %d, p = %d, C[N*N-1]=%f, wctime = %.4f\n",rank, N, procNum, C[N*N-1],wctime);
        }
        if(N <= 18){
            printMat(N,&C);
        }
    }
    //if im worker for the master
    else{
        barrier(); //wait for everyone to be ready before starting        
        wctime = 0.0;
        sizes[0]=MAT_SIZE;
        p[0]=NUM_PROCESSORS;   
         
        for (run=0; run<1; run++) {
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
            
            //initialize A and C row block
            sizeA = calcSize(rank, blockSize);
            initRowBlk(sizeA, &ArowBlock);
            sizeB = sizeA;
            initColBlk(sizeB, &BcolBlock);
            sizeC=blockSize*N;
            initRowBlk(sizeC, &CrowBlock);
            //receive my own A row and B column block 
            MPI_Irecv(BcolBlock, sizeB, MPI_DOUBLE, 0, rank+0xf000, MPI_COMM_WORLD, &request);
            MPI_Irecv(ArowBlock, sizeA, MPI_DOUBLE, 0, rank+0xe000, MPI_COMM_WORLD, &requestBackup);
            
            //must wait for receive to complete,no shortcut here 
            MPI_Wait(&request, &status);
            MPI_Wait(&requestBackup, &status);
            //self calculation
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock, &BcolBlock, &CrowBlock);


            srcRank=rank-1;//this is for receiving B blocks
            sendTag=rank;//this is for sending B block              
            //we need to wrap around the message passing
            destRank = (rank == procNum -1) ? 0 : rank+1;
            while(srcRank!=rank){
                
                sizeT = calcSize(srcRank, blockSize); 
                BtempBlock = (double *)calloc(sizeT,sizeof(double));
                
                //recall: I use rank # as my tag 
                MPI_Irecv(BtempBlock, sizeT, MPI_DOUBLE,rank-1,srcRank,MPI_COMM_WORLD,&request);
                
                MPI_Wait(&request, &status);
                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                //we dont send this B back to the original owner
                if(destRank != sendTag){
                    MPI_Isend(BcolBlock, sizeB, MPI_DOUBLE, destRank, sendTag, MPI_COMM_WORLD,&request);
                }
                
                //do calc based on newly arrived block
                //parallelized matmul and send for larger blocks
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                //wait for send to complete then we can free the buffer 
                MPI_Wait(&request, &status);
                
                if(destRank != sendTag){
                    free(BcolBlock);
                    
                    sizeB = sizeT;
                    BcolBlock = (double *)calloc(sizeB,sizeof(double));
                    memcpy(BcolBlock, BtempBlock, sizeB*sizeof(double));
                    updateIndx(&srcRank,procNum);
                    updateIndx(&sendTag,procNum);                
                    
                    free(BtempBlock);
                }
                
            }
            
            //collective MPI    
            C = (double *)malloc(N*N*sizeof(double));
            gather(&CrowBlock, sizeC, &C); 
            printf("proc %d: N = %d, p = %d, wctime = %.4f\n",rank, N, procNum, wctime);
        }
    }
    
    MPI_Finalize();           

}
