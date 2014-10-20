#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int i, wctime;
   int ts, thisbody, otherbody;
   double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
   double etime, etime0, etime1, cptime;
   Body *oct, *wildCardsTo, *wildCardsFrom;
   Body myOct;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
        //read data from file
        //raw data is defined in this function
        initData();

        //catagolize bodies into octants
        //oct is defined in this function
        sliceOctants(&oct);
        
        //wait for everyone to be ready before starting timer
        //then boardcast N, K, Dt
        MPI_Barrier(MPI_COMM_WORLD); 
        boardcastConsts();
        
        //prepare buffers to scat octants' bodies
        prepScat(&oct);

        //distribute Octants to their owners
        //myOct is defined in this function
        MPI_Barrier(MPI_COMM_WORLD); 
        scatOctants(&myOct);
        
        //free Octants Buffer after scatting the boddies, 
        //free initial data
        freeOctants(&oct);
        freeInitData();

        //debug prints
        printOct(&myOct,rank);

        //estimate DU of each of my bodies
        //This defines what wildCardsTo is (which wildcard to send to dest).
        estimateDU(&myOct, &wildCardsTo);

        exchangeCards(&wildCardsTo, &wildCardsFrom);
   }else{
        //wait for everyone to be ready before starting timer
        //waiting for master's boardcast
        MPI_Barrier(MPI_COMM_WORLD); 
        boardcastConsts();
        
        //distribute Octants to their owners
        MPI_Barrier(MPI_COMM_WORLD); 
        scatOctants(&myOct);

        //debug prints
        //printOct(&myOct,rank);
        
        //estimate DU of each of my bodies
        estimateDU(&myOct, &wildCardsTo);
        
        exchangeCards(&wildCardsTo, &wildCardsFrom);
   }
   MPI_Finalize();           
}
