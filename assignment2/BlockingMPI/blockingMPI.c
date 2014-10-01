#include "mpi.h"
#include "matmul.h"

#define MAT_SIZE 80
#define NUM_PROCESSORS 4

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
            N = sizes[run];//matrix size
            blockSize = N/p[run];//# of rows/cols per block A/:q
           
            //init A row block and C row block for this process (rank 0)
            //C row block has constant size, calculated insize initRowBlk()
            sizeA = calcSize(rank, blockSize);
            sizeC = blockSize * N; 
            
            //init row block A and C based on sizeA and sizeC
            initRowBlk(sizeA, sizeC, &ArowBlock, &CrowBlock);
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
            
            //this B Column block belongs to master, so the tag is master's rank
            sizeB = calcSize(rank, blockSize); 
            initColBlk(sizeB, &BcolBlock);
            
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
            printf("%d finished self calc\n",rank);            
            
            MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
            printf("%d sending a B to %d with tag %d\n",rank, rank+1, rank); 
            free(BcolBlock); //up to here, B column block buffer has done its job
            
            for(srcRank=procNum-1; srcRank>0; srcRank--){                 
                printf("%d waiting for  %d's B\n",rank, srcRank);
                sizeT = calcSize(srcRank, blockSize); 
                BtempBlock = (double *)malloc(sizeT*sizeof(double));
                printf("%d predicts sizeB = %d\n",rank, sizeT);
                //recall: I use rank # as my tag, but master only receives message from last worker
                MPI_Recv(BtempBlock, sizeT, MPI_DOUBLE,procNum-1,srcRank,MPI_COMM_WORLD,&status);
                printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
                printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
                printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
                
                
                //NOTE the differences here, im passing Temp B buffer instead of B Column buff
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                //we dont send this B back to the original owner
                if(srcRank != rank+1){
                    printf("%d sending a B to %d with tag %d\n",rank, rank+1, srcRank); 
                    MPI_Send(BtempBlock, sizeT, MPI_DOUBLE, rank+1, srcRank, MPI_COMM_WORLD);
                }

                printf("d");
                free(BtempBlock);
                printf("e\n");
            }
            printf("%d completed\n",rank);
            
            C = (double *)malloc(N*N*sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
            MPI_Gather(CrowBlock,sizeC, MPI_DOUBLE, C, sizeC, MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
            //printf("%d spent %.2f seconds to get val = %.2f", rank, wctime, C[N*N-1]);
        }
        printMat(N,&C);
    }
    //if im worker for the master
    else{
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
        
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
            
            //variable sizeAB due to different row and col blocks each processor handles
            sizeB = calcSize(rank, blockSize); 
            BcolBlock = (double *)malloc(sizeB*sizeof(double));
            
            //receive my own B column block 
            MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);
            printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
            printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
            printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
            printf("%d received its own B\n",rank);           
            
            //self calculation
            wctime += matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock, &BcolBlock, &CrowBlock);
            printf("%d finished self calc\n",rank);            

            srcRank=rank-1;
            sendTag=rank;              
            while(srcRank!=rank){
                
                printf("%d waiting for  %d's B\n",rank, srcRank);
                sizeT = calcSize(srcRank, blockSize); 
                BtempBlock = (double *)malloc(sizeT*sizeof(double));
                
                printf("%d predicts sizeB = %d\n",rank, sizeT);
                //recall: I use rank # as my tag 
                MPI_Recv(BtempBlock, sizeT, MPI_DOUBLE,rank-1,srcRank,MPI_COMM_WORLD,&status);
                printf("%d status.Error = %d\n",rank,status.MPI_ERROR);
                printf("%d status.MPI_TAG = %d\n",rank,status.MPI_TAG);
                printf("%d status.MPI_SOURCE = %d\n",rank,status.MPI_SOURCE);
                
                //do calc based on newly arrived block
                wctime += matmul(srcRank,N,blockSize,sizeA,sizeT,&ArowBlock,&BtempBlock,&CrowBlock);
                //we need to wrap around the message passing
                destRank = (rank == procNum -1) ? 0 : rank+1;
                
                if(destRank != sendTag){
                    //this B Column block belongs to MPI src, so the tag is MPI src's rank
                    //we dont send this B back to the original owner
                    printf("%d sending a B to %d with tag %d\n",rank, destRank, sendTag); 
                    MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, destRank, sendTag, MPI_COMM_WORLD); 
                    free(BcolBlock);
                }
                
                sizeB = sizeT;
                BcolBlock = (double *)malloc(sizeB*sizeof(double));
                memcpy(BcolBlock, BtempBlock, sizeB);
                
                if(srcRank != 0){
                    srcRank--;
                }else{
                    srcRank = procNum-1;
                }
                
                if(sendTag != 0){
                    sendTag--;
                }else{
                    sendTag = procNum-1;
                }
                
                free(BtempBlock);
            }
            printf("%d completed\n",rank);
            
            C = (double *)malloc(N*N*sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
            MPI_Gather(CrowBlock,sizeC, MPI_DOUBLE, C, sizeC, MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
        }
    }
    
    MPI_Finalize();           

}
