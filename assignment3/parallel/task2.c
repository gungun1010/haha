#include "calc.h"
#include "mpi.h"

int i, wctime;
int ts;
Body *oct, *wildCardsTo;
Body myOct, myWildCards, newComer;
double etime, etime0, etime1, cptime;

//top-level function to wrap around body movement operations, called in main()
void timeLapse(){
    //if(rank == ROOT)printf("**************************** %d ***************************\n",ts);
    barrier();                                                     
    if (ts%128 == 0){
         timing(&etime0,&cptime);
         output(ts, &myOct); // Print output if necessary
         timing(&etime1, &cptime);
         etime = etime + (etime1 - etime0);
    }
    //estimate DU of each of my bodies
    //This defines what wildCardsTo is (which wildcard to send to dest).
    barrier();

    timing(&etime0,&cptime);
    estimateDU(&myOct, &wildCardsTo);
    timing(&etime1, &cptime);
    etime = etime + (etime1 - etime0);
    //scat my wildcards to all 
    barrier();
    prepScatWildcards(&wildCardsTo);
    barrier();
    exchangeCards(&myWildCards);
    barrier();
    //main force calculation 
    timing(&etime0,&cptime);
    calcForce(&myOct, &myWildCards);
    timing(&etime1, &cptime);
    etime = etime + (etime1 - etime0);
    
    //must free all buffers used in collective operations
    barrier();
    freeBuffer(WILDCARD_STAGE);
    barrier();
    //check if any bodies is leaving my octant
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

main(int argc, char **argv){
   double etime, etime0, etime1, cptime;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
        etime=0;
        /////////////////////////////////////////////////////////
        //                    init  Phase                      //
        /////////////////////////////////////////////////////////
        
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
        freeBuffer(INIT_STAGE);
        //debug prints
        //printOct(&myOct,rank);
        
        /////////////////////////////////////////////////////////
        //                    Time Step Phase                  //
        /////////////////////////////////////////////////////////
        for(ts=0; ts<K; ts++){
            timeLapse();
        }

        timing(&etime0,&cptime);
        output(ts, &myOct); // Print output if necessary
        timing(&etime1, &cptime);
        etime = etime + (etime1 - etime0);
        
        barrier(); 
        printf ("\nProc %d: Time for %d timesteps with %d bodies: %.3le seconds\n",rank, K, N, etime);
   }else{
        etime=0;
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
            timeLapse();
        }

        timing(&etime0,&cptime);
        output(ts, &myOct); // Print output if necessary
        timing(&etime1, &cptime);
        etime = etime + (etime1 - etime0);
        
        barrier();
        printf ("\nProc %d: Time for %d timesteps with %d bodies: %.3le seconds\n",rank, K, N, etime);
   }
   MPI_Finalize();           
}
