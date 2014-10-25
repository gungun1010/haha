#include "calc.h"
//#include "comm.h"
//#include "functions.h"
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
        //printOct(&myOct,rank);
        
        /////////////////////////////////////////////////////////
        //                    Time Step Phase                  //
        /////////////////////////////////////////////////////////
        for(ts=0; ts<K; ts++){
            printf("**************************** %d ***************************\n",ts);
            barrier();                                                     
            if (ts%128 == 0) output(ts, &myOct); // Print output if necessary
            //estimate DU of each of my bodies
            //This defines what wildCardsTo is (which wildcard to send to dest).
            //comp. func.
            barrier();
            estimateDU(&myOct, &wildCardsTo);

            //scat my wildcards to all 
            //comm. func.
            barrier();
            prepScatWildcards(&wildCardsTo);
            barrier();
            exchangeCards(&myWildCards);
            barrier();
            //comp. func.
            calcForce(&myOct, &myWildCards);
            
            //printOct(&myOct,rank);
            //must free all buffers used in collective operations
            barrier();
            freeBuffer(WILDCARD_STAGE);
            //comp. func.
            barrier();
            updateOwner(&oct, &myOct);//recycle the oct buffer here
            
            //prepare buffer for newcomers 
            barrier();
            prepScatNewcomer(&oct);
            //exchange newcomers for each octant
            barrier();
            exchangeNewcomer(&newComer);

            //append new comers into my octant
            barrier();
            welcomeNewcomer(&myOct, &newComer);
            //free all kinds of buffers
            barrier();
            freeCardsDeck(&wildCardsTo);
            barrier();
            freeNewcomer(&newComer);
            barrier();
            freeWildCards(&myWildCards);
            barrier();
            freeOctants(&oct);
            //must free all buffers used in collective operations
            barrier();
            freeBuffer(NEWCOMER_STAGE);
        }
        output(ts, &myOct); // Print output if necessary
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
        //printOct(&myOct,rank);
        
        /////////////////////////////////////////////////////////
        //                    Time Step Phase                  //
        /////////////////////////////////////////////////////////
         
        for(ts=0; ts<K; ts++){
            barrier();  
            if (ts%128 == 0) output(ts, &myOct); // Print output if necessary
            //estimate DU of each of my bodies
            barrier();
            estimateDU(&myOct, &wildCardsTo);

            //scat my wildcards to all 
            printf("2. exchange wildcards at %d\n",rank);
            barrier();
            prepScatWildcards(&wildCardsTo);
            barrier();
            exchangeCards(&myWildCards);
            
            barrier();
            calcForce(&myOct, &myWildCards);

            //printOct(&myOct,rank);
            //must free all buffers used in collective operations
            barrier();
            freeBuffer(WILDCARD_STAGE);
            barrier();
            updateOwner(&oct, &myOct);
            
            //prepare buffer for newcomers 
            barrier();
            prepScatNewcomer(&oct);

            //exchange newcomers for each octant
            barrier();
            exchangeNewcomer(&newComer);
            //FIXME: when there is a new comer program hangs
            //append new comers into my octant
            barrier();
            welcomeNewcomer(&myOct, &newComer);
            barrier();
            //free all kinds of buffers 
            freeCardsDeck(&wildCardsTo);
            barrier();
            freeNewcomer(&newComer);
            barrier();
            freeWildCards(&myWildCards);
            barrier();
            freeOctants(&oct);
            //must free all buffers used in collective operations
            barrier();
            freeBuffer(NEWCOMER_STAGE);
        }
        output(ts, &myOct); // Print output if necessary
   }
   MPI_Finalize();           
}
