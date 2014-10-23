#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int i, wctime;
   int ts;
   double etime, etime0, etime1, cptime;
   Body *oct, *wildCardsTo;
   Body myOct, myWildCards;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
        //read data from file
        //raw data is defined in this function
        //computation function
        initData();

        //catagolize bodies into octants
        //oct is defined in this function
        //comm. func.
        sliceOctants(&oct);
        
        //wait for everyone to be ready before starting timer
        //then boardcast N, K, Dt
        //comm. func.
        barrier(); 
        boardcastConsts();
        
        //prepare buffers to scat octants' bodies
        //comm. func.
        prepScat(&oct);

        //distribute Octants to their owners
        //myOct is defined in this function
        //comm. func.
        barrier(); 
        scatOctants(&myOct);
        
        //free Octants Buffer after scatting the boddies, 
        //free initial data
        //comm. func.
        freeOctants(&oct);
        //comp. func.
        freeInitData();
        //comm. func.
        freeBuffer();
        //debug prints
        //comp. func.
        //printOct(&myOct,rank);

        //estimate DU of each of my bodies
        //This defines what wildCardsTo is (which wildcard to send to dest).
        //comp. func.
        estimateDU(&myOct, &wildCardsTo);

        //scat my wildcards to all 
        //comm. func.
        prepScatWildcards(&wildCardsTo);
        exchangeCards(&myWildCards);
        
        //must free all buffers used in collective operations
        barrier(&oct, &myOct);
        //comp. func.
        calcForce(&myOct, &myWildCards);
        
        printOct(&myOct,rank);
        //comp. func.
        updateOwner(&oct, &myOct);//recycle the oct buffer here
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
        
        //must free all buffers used in collective operations
        barrier();
        calcForce(&myOct, &myWildCards);

        printOct(&myOct,rank);
        updateOwner(&oct, &myOct);
   }
   MPI_Finalize();           
}
