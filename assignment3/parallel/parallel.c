#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int i, wctime;
   int ts, thisbody, otherbody;
   double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
   double etime, etime0, etime1, cptime;
   Body *oct;
   Body myOct, wildCard;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
        //read data from file
        initData();
        //catagolize bodies into octants 
        sliceOctants(&oct);
        printf("finish init data\n"); 
        
        //wait for everyone to be ready before starting timer
        //then boardcast N, K, Dt
        MPI_Barrier(MPI_COMM_WORLD); 
        boardcastConsts();
        //we dont need init Size array after boardcast, free it free(initOctSize); 
        free(initOctSize);  
        
        //distribute Octants to their owners
        scatOctants(&oct);
        
        //receive the bodies within my octant
        recvBodies(&myOct);
        
        //free Octants Buffer after scatting the boddies, 
        //free initial data
        freeOctants(&oct);
        freeInitData();

        //debug prints
        //printOct(&myOct,rank);

        //estimate DU of each of my bodies
        estimateDU(&myOct, &wildCard);
   }else{
        //wait for everyone to be ready before starting timer
        //waiting for master's boardcast
        MPI_Barrier(MPI_COMM_WORLD); 
        boardcastConsts();
        
        //receive the bodies within my octant
        recvBodies(&myOct);

        //debug prints
        //printOct(&myOct,rank);
        
        //estimate DU of each of my bodies
        //esitmateDU(&myOct);
   }
   MPI_Finalize();           
}
