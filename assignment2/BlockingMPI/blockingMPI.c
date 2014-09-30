#include "mpi.h"
#include "matmul.h"

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
        sizes[0]=MAT_SIZE;
        p[0]=NUM_PROCESSORS;   
        
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting timer
        
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
            
            //variable sizeAB due to different row and col blocks each processor handles
            //after MPI_Send, work on my own task, rank=0
            sizeB = calcSize(rank, blockSize); 
            initColBlk(sizeB, &BcolBlock);

            wctime = matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
            
            //this B Column block belongs to master, so the tag is master's rank
            MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
            free(BcolBlock);
            
             
            for(srcRank=procNum-1; srcRank>0; srcRank--){                 
                sizeB = calcSize(srcRank, blockSize); 
                initColBlk(sizeB, &BcolBlock);
                
                //recall: I use rank # as my tag, but master only receives message from last worker
                MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE,procNum-1,srcRank,MPI_COMM_WORLD,&status);
                
                wctime = matmul(srcRank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, rank+1, srcRank, MPI_COMM_WORLD);
                free(BcolBlock);
            }


            //FIXME debug prints
            printf("rank %d :",rank);
            for(i=0; i<sizeC; i++)
                printf(" %.2f ",CrowBlock[i]);
            printf("\n");
        }
    }
    //if im worker for the master
    else{
        MPI_Barrier(MPI_COMM_WORLD); //wait for everyone to be ready before starting
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
            BcolBlock = (double *)calloc(sizeB, sizeof(double));
            
            //recall: I use rank # as my tag 
            MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);
                       
            wctime = matmul(rank,N,blockSize,sizeA,sizeB,&ArowBlock, &BcolBlock, &CrowBlock);
            
            //we need to wrap around the message passing
            destRank = rank == procNum -1 ? 0 : rank+1;
            
            //this B Column block belongs to this worker, so the tag is worker's rank
            MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, destRank, rank, MPI_COMM_WORLD); 
            free(BcolBlock);
            
            srcRank=rank-1;              
            while(srcRank!=rank){    
                sizeB = calcSize(srcRank, blockSize); 
                initColBlk(sizeB, &BcolBlock);
                
                //recall: I use rank # as my tag 
                MPI_Recv(BcolBlock, sizeB, MPI_DOUBLE,rank-1,srcRank,MPI_COMM_WORLD,&status);

                wctime = matmul(srcRank,N,blockSize,sizeA,sizeB,&ArowBlock,&BcolBlock,&CrowBlock);
                //we need to wrap around the message passing
                destRank = rank == procNum -1 ? 0 : rank+1;

                //this B Column block belongs to MPI src, so the tag is MPI src's rank
                MPI_Send(BcolBlock, sizeB, MPI_DOUBLE, destRank, srcRank, MPI_COMM_WORLD);
                free(BcolBlock);
                if(srcRank != 0){
                    srcRank--;
                }else{
                    srcRank = procNum-1;
                }
            }
            //FIXME debug prints
            printf("rank %d :",rank);
            for(i=0; i<sizeC; i++)
                printf(" %.2f ",CrowBlock[i]);
            printf("\n");
            
        }
    }
    
    MPI_Finalize();           

}
