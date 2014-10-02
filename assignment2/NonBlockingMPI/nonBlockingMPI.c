#include "mpi.h"
#include "matmul.h"
#include "mpiWrappers.h"

#define MAT_SIZE 8
#define NUM_PROCESSORS 8

main(int argc, char **argv) {
    int N, i, run, blockSize, sizeA,sizeB,sizeC,sizeT;
    double *C,*ArowBlock,*BtempBlock,*BcolBlock,*CrowBlock;
    int sizes[1];//matrix size
    int p[1];//number of processors 
    double wctime, sparm;
    int rank, procNum;
    int srcRank, destRank,sendTag;
    double worktime;
    MPI_Status status;  
    MPI_Request request;
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
            
            //initialize A and C row block
            initAnC(rank, blockSize, N, &sizeA, &sizeC, &ArowBlock, &CrowBlock);

            //distribute B Col block to each of the workers
            distributeB(rank, procNum, blockSize, &BcolBlock,&request);
            
            //this B Column block belongs to master, so the tag is master's rank
            sizeB = calcSize(rank, blockSize); 
            initColBlk(sizeB, &BcolBlock);
           
            MPI_Isend(BcolBlock, sizeB, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD, &request);
            //after MPI_Isend, work on my own task, rank=0
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
            
            free(BcolBlock); //up to here, B column block buffer has done its job
            MPI_Wait(&request, &status);
             
            for(srcRank=procNum-1; srcRank>0; srcRank--){                 
                sizeT = calcSize(srcRank, blockSize); 
                BtempBlock = (double *)calloc(sizeT,sizeof(double));
                
                //recall: I use rank # as my tag, but master only receives message from last worker
                MPI_Irecv(BtempBlock, sizeT, MPI_DOUBLE,procNum-1,srcRank,MPI_COMM_WORLD,&request);
                MPI_Wait(&request, &status);
                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                //we dont send this B back to the original owner
                if(srcRank != rank+1){
                    MPI_Isend(BtempBlock, sizeT, MPI_DOUBLE, rank+1, srcRank, MPI_COMM_WORLD, &request);
                }//NOTE the differences here, im passing Temp B buffer instead of B Column buff
                
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                free(BtempBlock);
                MPI_Wait(&request, &status);
            }
            
            //collective MPI 
            C = (double *)malloc(N*N*sizeof(double));
            gather(&CrowBlock, sizeC, &C); 
            printf("proc %d: N = %d, p = %d, C[N*N-1]=%f, wctime = %.4f\n",rank, N, procNum, C[N*N-1],wctime);
        }
        printMat(N,&C);
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
            initAnC(rank, blockSize, N, &sizeA, &sizeC, &ArowBlock, &CrowBlock);
            
            //variable sizeAB due to different row and col blocks each processor handles
            sizeB = calcSize(rank, blockSize); 
            BcolBlock = (double *)calloc(sizeB,sizeof(double));
            
            //receive my own B column block 
            MPI_Irecv(BcolBlock, sizeB, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request);
            
            MPI_Wait(&request, &status);
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
                MPI_Isend(BcolBlock, sizeB, MPI_DOUBLE, destRank, sendTag, MPI_COMM_WORLD,&request);
                
                //do calc based on newly arrived block
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                
                if(destRank != sendTag){
                    free(BcolBlock);
                    
                    sizeB = sizeT;
                    BcolBlock = (double *)calloc(sizeB,sizeof(double));
                    memcpy(BcolBlock, BtempBlock, sizeB*sizeof(double));
                    updateIndx(&srcRank,procNum);
                    updateIndx(&sendTag,procNum);                
                    
                    free(BtempBlock);
                }
                
                MPI_Wait(&request, &status);
            }
            
            //collective MPI    
            C = (double *)malloc(N*N*sizeof(double));
            gather(&CrowBlock, sizeC, &C); 
            printf("proc %d: N = %d, p = %d, wctime = %.4f\n",rank, N, procNum, wctime);
        }
    }
    
    MPI_Finalize();           

}
