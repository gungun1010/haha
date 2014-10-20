#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int i, wctime;
   int ts, thisbody, otherbody;
   double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
   double etime, etime0, etime1, cptime;
   Body *oct, *wildCardsTo;
   Body myOct, myWildCards;
    
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
        barrier(); 
        boardcastConsts();
        
        //prepare buffers to scat octants' bodies
        prepScat(&oct);

        //distribute Octants to their owners
        //myOct is defined in this function
        barrier(); 
        scatOctants(&myOct);
        
        //free Octants Buffer after scatting the boddies, 
        //free initial data
        freeOctants(&oct);
        freeInitData();
        freeBuffer();
        //debug prints
        printOct(&myOct,rank);

        //estimate DU of each of my bodies
        //This defines what wildCardsTo is (which wildcard to send to dest).
        estimateDU(&myOct, &wildCardsTo);

        //scat my wildcards to all 
        prepScatWildcards(&wildCardsTo);
        exchangeCards(&myWildCards);
   }else{
        //wait for everyone to be ready before starting timer
        //waiting for master's boardcast
        barrier(); 
        boardcastConsts();
        
        //distribute Octants to their owners
        barrier(); 
        scatOctants(&myOct);

        //debug prints
        //printOct(&myOct,rank);
        
        //estimate DU of each of my bodies
        estimateDU(&myOct, &wildCardsTo);

        //scat my wildcards to all 
        prepScatWildcards(&wildCardsTo);
        exchangeCards(&myWildCards);
   }
   MPI_Finalize();           
}
