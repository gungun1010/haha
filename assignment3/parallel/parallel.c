#include "functions.h"
#include "mpi.h"

main(int argc, char **argv){
   int wctime;
   int ts, thisbody, otherbody;
   double vavgx, vavgy, vavgz, ax, ay, az, deltaf[3];
   double etime, etime0, etime1, cptime;
    
   init(&argc, &argv, &procNum, &rank);
   
   if (rank == 0){
      
        //read data from file
        initData();
        printf("finish init data\n"); 
        //wait for everyone to be ready before starting timer
        MPI_Barrier(MPI_COMM_WORLD); 
        //slice them into octants
        boardcastConsts();
        
        //catagolize bodies into octants 
        scatOctants();
   
   }else{
        //wait for everyone to be ready before starting timer
        MPI_Barrier(MPI_COMM_WORLD); 
        //slice them into octants
        boardcastConsts();
        
   }
   MPI_Finalize();           
}
