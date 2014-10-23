#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int i, wctime;
   int ts;
   double etime, etime0, etime1, cptime;
   Body *oct, *wildCardsTo;
   Body myOct, myWildCards, newComer;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
        /////////////////////////////////////////////////////////
        //                    init  Phase                      //
        /////////////////////////////////////////////////////////
        
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
        freeBuffer(INIT_STAGE);
        //debug prints
        //comp. func.
        printOct(&myOct,rank);
        
        /////////////////////////////////////////////////////////
        //                    Time Step Phase                  //
        /////////////////////////////////////////////////////////
                                                       
        //estimate DU of each of my bodies
        //This defines what wildCardsTo is (which wildcard to send to dest).
        //comp. func.
        estimateDU(&myOct, &wildCardsTo);

        //scat my wildcards to all 
        //comm. func.
        prepScatWildcards(&wildCardsTo);
        exchangeCards(&myWildCards);
        //FIXME: need to free wildcard 
        barrier(&oct, &myOct);
        //comp. func.
        calcForce(&myOct, &myWildCards);
        
        //printOct(&myOct,rank);
        //must free all buffers used in collective operations
        freeBuffer(WILDCARD_STAGE);
        //comp. func.
        updateOwner(&oct, &myOct);//recycle the oct buffer here
        
        //prepare buffer for newcomers 
        prepScatNewcomer(&oct);

        //FIXME: implement this and then we are there, hopefully 
        //exchange newcomers for each octant
        exchangeNewcomer(&newComer);

        //append new comers into my octant
        welcomeNewcomer(&myOct, &newComer);
   }else{
        /////////////////////////////////////////////////////////
        //                    init  Phase                      //
        /////////////////////////////////////////////////////////

        //wait for everyone to be ready before starting timer
        //waiting for master's boardcast
        barrier(); 
        boardcastConsts();
        
        //distribute Octants to their owners
        barrier(); 
        scatOctants(&myOct);

        //debug prints
        printOct(&myOct,rank);
        
        /////////////////////////////////////////////////////////
        //                    Time Step Phase                  //
        /////////////////////////////////////////////////////////
        
        //estimate DU of each of my bodies
        estimateDU(&myOct, &wildCardsTo);

        //scat my wildcards to all 
        prepScatWildcards(&wildCardsTo);
        exchangeCards(&myWildCards);
        
        barrier();
        calcForce(&myOct, &myWildCards);

        //printOct(&myOct,rank);
        //must free all buffers used in collective operations
        freeBuffer(WILDCARD_STAGE);
        updateOwner(&oct, &myOct);
        
        //prepare buffer for newcomers 
        prepScatNewcomer(&oct);

        //exchange newcomers for each octant
        exchangeNewcomer(&newComer);

        //append new comers into my octant
        welcomeNewcomer(&myOct, &newComer);
   }
   MPI_Finalize();           
}
