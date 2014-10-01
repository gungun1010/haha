#include "mpi.h"
#include "matmul.h"

#define MAT_SIZE 68
#define NUM_PROCESSORS 2

main(int argc, char **argv) {
    int N, i, run, blockSize, sizeA,sizeB,sizeC;
    double *A, *B, *C, *ArowBlock, *BcolBlock, *CrowBlock;

    int sizes[1];//matrix size
    int p[1];//number of processors 
    double wctime, sparm;
    int rank, procNum;
    int srcRank, destRank;
    double worktime;
    MPI_Status status;  

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&procNum); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

    /* If I am the master (rank 0) ... */
    if (rank == 0) {
        //sparm = 0; //initialize the workers' work times 
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
        wctime = 0.0;
        sizes[0]=MAT_SIZE;
        p[0]=NUM_PROCESSORS;   
        
        
        for (run=0; run<1; run++) {
            printf("%d start\n",rank);     
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows/cols per block A/:q
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * N; 
            
            //init row block A and C based on sizeA and sizeC
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            printf("%d finished initRowBlk\n",rank);
            //send B column block (tag = i) to each process with rank i
            for(i=rank;i<procNum;i++){

                //variable sizeAB due to different row and col blocks each processor handles
                sizeB = calcSize(i, blockSize); 
                initColBlk(sizeB, &BcolBlock);
                
                //using target rank as tag
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
                
                //recycle B column block memories
                free(BcolBlock);
            }
            
            printf("%d finished sending block B\n",rank);
            //variable sizeAB due to different row and col blocks each processor handles
            //after MPI_Send, work on my own task, rank=0
            sizeB = calcSize(rank, blockSize); 
            initColBlk(sizeB, &BcolBlock);

            MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
            printf("%d finished sending self B(%d)\n",rank,sizeB);
            
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
            printf("%d finished self calc\n",rank);            
            //this B Column block belongs to master, so the tag is master's rank
            free(BcolBlock);

            
            for(srcRank=procNum-1; srcRank>0; srcRank--){                 
                printf("%d waiting for  %d's B\n",rank, srcRank);
                sizeB = calcSize(srcRank, blockSize); 
                BcolBlock = (double *)malloc(sizeB*sizeof(double));
                printf("%d predicts sizeB = %d\n",rank, sizeB);
                //recall: I use rank # as my tag, but master only receives message from last worker
                MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE,procNum-1,srcRank,MPI_COMM_WORLD,&status);
                printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
                printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
                printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
                printf("%d b ",rank);
                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                if(srcRank != rank+1){
                    MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, srcRank, MPI_COMM_WORLD);
                    printf("c");
                }

                wctime += matmul(srcRank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
                printf("d");
                free(BcolBlock);
                printf("e\n");
            }
            printf("%d completed\n",rank);
            
            C = (double *)malloc(N*N*sizeof(double));
            MPI_Gather(CrowBlock,sizeC, MPI_DOUBLE, C, sizeC, MPI_DOUBLE,0,MPI_COMM_WORLD);
            printf("%d spent %.2f seconds to get val = %.2f", rank, wctime, C[N*N-1]);
        }
        //printMat(N,&C);
    }
    //if im worker for the master
    else{
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
        
        printf("%d start\n",rank);     
        wctime = 0.0;
        sizes[0]=MAT_SIZE;
        p[0]=NUM_PROCESSORS;   
        
        for (run=0; run<1; run++) {
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows and cols per block A and B
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * N; 
            
            //init row block A and C based on size
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
            
            printf("%d finished initRowBlk\n",rank);
            //variable sizeAB due to different row and col blocks each processor handles
            sizeB = calcSize(rank, blockSize); 
            BcolBlock = (double *)malloc(sizeB*sizeof(double));
            
            //recall: I use rank # as my tag 
            MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);
            printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
            printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
            printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
            printf("%d received its own B\n",rank);           
            //we need to wrap around the message passing
            destRank = (rank == procNum -1) ? 0 : rank+1;
            printf("%d sending own B to %d\n",rank, destRank); 
            //this B Column block belongs to this worker, so the tag is worker's rank
            MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, destRank, rank, MPI_COMM_WORLD); 
            printf("%d sending a B to %d with tag %d \n",rank,destRank,rank);
            
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock, &BcolBlock, &CrowBlock);
            
            printf("%d finished self calc\n",rank);            
            free(BcolBlock);
            
            srcRank=rank-1;              
            while(srcRank!=rank){    
                printf("%d waiting for  %d's B\n",rank, srcRank);
                sizeB = calcSize(srcRank, blockSize); 
                BcolBlock = (double *)malloc(sizeB*sizeof(double));
                
                printf("%d predicts sizeB = %d\n",rank, sizeB);
                //recall: I use rank # as my tag 
                MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE,rank-1,srcRank,MPI_COMM_WORLD,&status);
                printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
                printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
                printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
                
                //we need to wrap around the message passing
                destRank = (rank == procNum -1) ? 0 : rank+1;
                if(destRank != srcRank){
                    //this B Column block belongs to MPI src, so the tag is MPI src's rank
                    MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, destRank, srcRank, MPI_COMM_WORLD);
                    
                    printf("%d sending a B to %d with tag %d \n",rank,destRank,srcRank);
                }
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
                printf("%d d ",rank);
                free(BcolBlock);
                printf("%d e\n",rank);
                if(srcRank != 0){
                    srcRank--;
                }else{
                    srcRank = procNum-1;
                }
            }
            printf("%d completed\n",rank);
            
            C = (double *)malloc(N*N*sizeof(double));
            MPI_Gather(CrowBlock,sizeC, MPI_DOUBLE, C, sizeC, MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
    }
    
    MPI_Finalize();           

}
