#include "comm.h"
#include "mpi.h"

void barrier(){
    MPI_Barrier(MPI_COMM_WORLD);
}

void initBody(Body *a, size_t initialSize) {
      a->x = (double *)malloc(initialSize * sizeof(double));
      a->y = (double *)malloc(initialSize * sizeof(double));
      a->z = (double *)malloc(initialSize * sizeof(double));
      a->vx = (double *)malloc(initialSize * sizeof(double));
      a->vy = (double *)malloc(initialSize * sizeof(double));
      a->vz = (double *)malloc(initialSize * sizeof(double));
      a->mass = (double *)malloc(initialSize * sizeof(double));
      
      a->used = 0;
      a->size = initialSize;
}

void insertBody(Body *a, int i){
      if (a->used == a->size) {
          a->size++;
          
          a->x = (double *)realloc(a->x, a->size * sizeof(double));
          a->y = (double *)realloc(a->y, a->size * sizeof(double));
          a->z = (double *)realloc(a->z, a->size * sizeof(double));
          a->vx = (double *)realloc(a->vx, a->size * sizeof(double));
          a->vy = (double *)realloc(a->vy, a->size * sizeof(double));
          a->vz = (double *)realloc(a->vz, a->size * sizeof(double));
          a->mass = (double *)realloc(a->mass, a->size * sizeof(double));
      }
      
      a->x[a->used] = x[i];
      a->y[a->used] = y[i];
      a->z[a->used] = z[i];
      a->vx[a->used] = vx[i];
      a->vy[a->used] = vy[i];
      a->vz[a->used] = vz[i];
      a->mass[a->used] = mass[i];
      a->used++;
}

void removeBody(Body *a, int index){
   int i;
   int newSize;

   newSize = a->size -1;
   for(i=index; i< newSize; i++){
       a->x[i] = a->x[i+1];
       a->y[i] = a->y[i+1];
       a->z[i] = a->z[i+1];

       a->vx[i] = a->vx[i+1];
       a->vy[i] = a->vy[i+1];
       a->vz[i] = a->vz[i+1];

       a->fx[i] = a->fx[i+1];
       a->fy[i] = a->fy[i+1];
       a->fz[i] = a->fz[i+1];

       a->mass[i] = a->mass[i+1];
    }
    
    a->x = (double *)realloc(a->x, newSize * sizeof(double));
    a->y = (double *)realloc(a->y, newSize * sizeof(double));
    a->z = (double *)realloc(a->z, newSize * sizeof(double));
    
    a->vx = (double *)realloc(a->vx, newSize * sizeof(double));
    a->vy = (double *)realloc(a->vy, newSize * sizeof(double));
    a->vz = (double *)realloc(a->vz, newSize * sizeof(double));
    
    a->fx = (double *)realloc(a->fx, newSize * sizeof(double));
    a->fy = (double *)realloc(a->fy, newSize * sizeof(double));
    a->fz = (double *)realloc(a->fz, newSize * sizeof(double));
    
    a->mass = (double *)realloc(a->mass, newSize * sizeof(double));

    a->size = a->used = newSize;
}

//MPI initialization function
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef){
    
    MPI_Init(argcRef, argvRef);

    MPI_Comm_size(MPI_COMM_WORLD, procNumRef); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, rankRef); // Which process am I?
}

void boardcastConsts(){
        
    //boardcast N, K, Dt
    MPI_Bcast(&N, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    //scat sizeArr of each octant for the owner
    if(rank !=0)
        sizeArr = calloc(procNum, sizeof(int));

    MPI_Bcast(sizeArr, procNum, MPI_INT, ROOT, MPI_COMM_WORLD);
}

void prepScat(Body** oct){
   int i;
   int bodyNum=0;
   int offset=0;

    //get total number of bodies in the problem
   for (i=0; i<procNum; i++)
       bodyNum = bodyNum + sizeArr[i];
   
   //allocate memory to buffers 
   xArr =(double*) calloc(bodyNum, sizeof(double));      
   yArr =(double*) calloc(bodyNum, sizeof(double));      
   zArr =(double*) calloc(bodyNum, sizeof(double));      
   
   vxArr =(double*) calloc(bodyNum, sizeof(double));      
   vyArr =(double*) calloc(bodyNum, sizeof(double));      
   vzArr =(double*) calloc(bodyNum, sizeof(double));      

   massArr =(double*) calloc(bodyNum, sizeof(double));      
   
   displ =(int*) calloc(procNum, sizeof(int));
   //append things into their array buffer
   for(i=0; i<procNum;i++){
       displ[i] = offset;
       memcpy(xArr+offset, (*oct)[i].x, (*oct)[i].size*sizeof(double));
       memcpy(yArr+offset, (*oct)[i].y, (*oct)[i].size*sizeof(double));
       memcpy(zArr+offset, (*oct)[i].z, (*oct)[i].size*sizeof(double));

       memcpy(vxArr+offset, (*oct)[i].vx, (*oct)[i].size*sizeof(double));
       memcpy(vyArr+offset, (*oct)[i].vy, (*oct)[i].size*sizeof(double));
       memcpy(vzArr+offset, (*oct)[i].vz, (*oct)[i].size*sizeof(double));

       memcpy(massArr+offset, (*oct)[i].mass, (*oct)[i].size*sizeof(double));

       offset+=(*oct)[i].size;
   }
}

void scatOctants(Body* myOct){
   //allocate coordinates for bodies in the octant 
   mySize = sizeArr[rank];
   
   myOct->x = (double *) calloc(mySize, sizeof(double)); 
   myOct->y = (double *) calloc(mySize, sizeof(double)); 
   myOct->z = (double *) calloc(mySize, sizeof(double)); 
   //allocate velocity for bodies in the octant
   myOct->vx = (double *) calloc(mySize, sizeof(double)); 
   myOct->vy = (double *) calloc(mySize, sizeof(double)); 
   myOct->vz = (double *) calloc(mySize, sizeof(double)); 
    
   //allocate masses for bodies in the octant
   myOct->mass = (double *) calloc(mySize, sizeof(double)); 
   
   //scat array with specified length 
   MPI_Scatterv(xArr, sizeArr, displ, MPI_DOUBLE, myOct->x, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   MPI_Scatterv(yArr, sizeArr, displ, MPI_DOUBLE, myOct->y, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   MPI_Scatterv(zArr, sizeArr, displ, MPI_DOUBLE, myOct->z, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

   MPI_Scatterv(vxArr, sizeArr, displ, MPI_DOUBLE, myOct->vx, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   MPI_Scatterv(vyArr, sizeArr, displ, MPI_DOUBLE, myOct->vy, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   MPI_Scatterv(vzArr, sizeArr, displ, MPI_DOUBLE, myOct->vz, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   
   MPI_Scatterv(massArr, sizeArr, displ, MPI_DOUBLE, myOct->mass, mySize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
   
   //each octant has its own init size, mySize 
   myOct->size = myOct->used = mySize;
}

void freeBuffer(int stage){
    free(massArr);
    free(xArr);
    free(yArr);
    free(zArr);

    //wild card stage only scats positions and mass
    //init stage has no force values to exchange
    if(stage == INIT_STAGE){
        free(vxArr);
        free(vyArr);
        free(vzArr);
    }

    //only OWNER stage exchanges forces
    if(stage == NEWCOMER_STAGE){
        free(vxArr);
        free(vyArr);
        free(vzArr);
        free(fxArr);
        free(fyArr);
        free(fzArr);
    }

    free(displ);
    free(sizeArr);

    if(stage == WILDCARD_STAGE || stage == NEWCOMER_STAGE){
        free(transDispl);
        free(transSizeArr);
    }
}

void prepScatWildcards(Body** wildCardsTo){
   int i;
   int bodyNum=0;
   int offset=0;
   int size;

    //get total number of bodies in the deck
   for (i=0; i<procNum; i++){
       if(i!=rank){
           bodyNum = bodyNum + (*wildCardsTo)[i].size;
        }
    }

   //allocate memory to buffers 
   xArr =(double*) calloc(bodyNum, sizeof(double));      
   yArr =(double*) calloc(bodyNum, sizeof(double));      
   zArr =(double*) calloc(bodyNum, sizeof(double));      

   massArr =(double*) calloc(bodyNum, sizeof(double));      
   
   displ =(int*) calloc(procNum, sizeof(int));
   sizeArr =(int*) calloc(procNum, sizeof(int));
   
   //append things into their array buffer
   //the Arr buffers collects all the elements that are wildcards to octants
   for(i=0; i<procNum;i++){
       if(i!=rank){
           size = (*wildCardsTo)[i].size;
           sizeArr[i] = size;
           displ[i] = offset;
           memcpy(xArr+offset, (*wildCardsTo)[i].x, size*sizeof(double));
           memcpy(yArr+offset, (*wildCardsTo)[i].y, size*sizeof(double));
           memcpy(zArr+offset, (*wildCardsTo)[i].z, size*sizeof(double));

           memcpy(massArr+offset, (*wildCardsTo)[i].mass, size*sizeof(double));

           offset+=size;
       }
   }

   //transpose initOctSize and displ
   transSizeArr =(int*) calloc(procNum, sizeof(int));
   transDispl =(int*) calloc(procNum, sizeof(int));
   barrier();
   MPI_Alltoall(sizeArr, 1, MPI_INT, transSizeArr, 1, MPI_INT, MPI_COMM_WORLD);
   
   offset=0;
   for(i=0;i<procNum;i++){
        transDispl[i] = offset;
        offset+=transSizeArr[i];
   }
}

void exchangeCards(Body* myWildCards){
    int i;

   //allocate coordinates for bodies in the octant 
   mySize = 0;
   for(i=0;i<procNum;i++){
       mySize += transSizeArr[i];
   }
    
   myWildCards->x = (double *) calloc(mySize, sizeof(double)); 
   myWildCards->y = (double *) calloc(mySize, sizeof(double)); 
   myWildCards->z = (double *) calloc(mySize, sizeof(double)); 
    
   //allocate masses for bodies in the octant
   myWildCards->mass = (double *) calloc(mySize, sizeof(double)); 
   
   //scat array with specified length 
   //each octant scat its wildcards' buffers
   barrier();
   //scat coordinates first
   MPI_Alltoallv(xArr, sizeArr, displ, MPI_DOUBLE, myWildCards->x, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
   MPI_Alltoallv(yArr, sizeArr, displ, MPI_DOUBLE, myWildCards->y, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
   MPI_Alltoallv(zArr, sizeArr, displ, MPI_DOUBLE, myWildCards->z, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
   
   //now scat mass
   MPI_Alltoallv(massArr, sizeArr, displ, MPI_DOUBLE, myWildCards->mass, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
   
   //each octant has its own wildcard size, mySize 
   myWildCards->size = myWildCards->used = mySize;
   //printf("%d: %d\n",rank, myWildCards->size);
}

void prepScatNewcomer(Body** oct){
   int i;
   int bodyNum=0;
   int offset=0;
   int size;

    //get total number of bodies in the deck
   for (i=0; i<procNum; i++){
       if(i!=rank){
           bodyNum = bodyNum + (*oct)[i].size;
        }
    }
   
   //allocate memory to buffers 
   xArr =(double*) calloc(bodyNum, sizeof(double));      
   yArr =(double*) calloc(bodyNum, sizeof(double));      
   zArr =(double*) calloc(bodyNum, sizeof(double));      

   //allocate memory to buffers 
   vxArr =(double*) calloc(bodyNum, sizeof(double));      
   vyArr =(double*) calloc(bodyNum, sizeof(double));      
   vzArr =(double*) calloc(bodyNum, sizeof(double));      
   
   //allocate memory to buffers 
   fxArr =(double*) calloc(bodyNum, sizeof(double));      
   fyArr =(double*) calloc(bodyNum, sizeof(double));      
   fzArr =(double*) calloc(bodyNum, sizeof(double));      
   
   massArr =(double*) calloc(bodyNum, sizeof(double));      
   
   displ =(int*) calloc(procNum, sizeof(int));
   sizeArr =(int*) calloc(procNum, sizeof(int));
   
   //append things into their array buffer
   //the Arr buffers collects all the elements that are wildcards to octants
   for(i=0; i<procNum;i++){
       if(i!=rank){
           size = (*oct)[i].size;
           sizeArr[i] = size;
           displ[i] = offset;
           memcpy(xArr+offset, (*oct)[i].x, size*sizeof(double));
           memcpy(yArr+offset, (*oct)[i].y, size*sizeof(double));
           memcpy(zArr+offset, (*oct)[i].z, size*sizeof(double));

           memcpy(vxArr+offset, (*oct)[i].vx, size*sizeof(double));
           memcpy(vyArr+offset, (*oct)[i].vy, size*sizeof(double));
           memcpy(vzArr+offset, (*oct)[i].vz, size*sizeof(double));
           
           memcpy(fxArr+offset, (*oct)[i].fx, size*sizeof(double));
           memcpy(fyArr+offset, (*oct)[i].fy, size*sizeof(double));
           memcpy(fzArr+offset, (*oct)[i].fz, size*sizeof(double));
           
           memcpy(massArr+offset, (*oct)[i].mass, size*sizeof(double));

           offset+=size;
       }
   }
   
   //transpose initOctSize and displ
   transSizeArr =(int*) calloc(procNum, sizeof(int));
   transDispl =(int*) calloc(procNum, sizeof(int));
   barrier();
   MPI_Alltoall(sizeArr, 1, MPI_INT, transSizeArr, 1, MPI_INT, MPI_COMM_WORLD);
   
   offset=0;
   for(i=0;i<procNum;i++){
        transDispl[i] = offset;
        offset+=transSizeArr[i];
   }
}

void exchangeNewcomer(Body* newComer){
    int i;
    int ncSize=0;
   
   //allocate coordinates for bodies in the octant
   for(i=0;i<procNum;i++){
       ncSize += transSizeArr[i];
   }
     
    newComer->x = (double *) calloc(ncSize, sizeof(double));
    newComer->y = (double *) calloc(ncSize, sizeof(double));
    newComer->z = (double *) calloc(ncSize, sizeof(double));

    newComer->vx = (double *) calloc(ncSize, sizeof(double));
    newComer->vy = (double *) calloc(ncSize, sizeof(double));
    newComer->vz = (double *) calloc(ncSize, sizeof(double));
    
    newComer->fx = (double *) calloc(ncSize, sizeof(double));
    newComer->fy = (double *) calloc(ncSize, sizeof(double));
    newComer->fz = (double *) calloc(ncSize, sizeof(double));
    
    newComer->mass = (double *) calloc(ncSize, sizeof(double));

    barrier();
    
    MPI_Alltoallv(xArr, sizeArr, displ, MPI_DOUBLE, newComer->x, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(yArr, sizeArr, displ, MPI_DOUBLE, newComer->y, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(zArr, sizeArr, displ, MPI_DOUBLE, newComer->z, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Alltoallv(vxArr, sizeArr, displ, MPI_DOUBLE, newComer->vx, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(vyArr, sizeArr, displ, MPI_DOUBLE, newComer->vy, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(vzArr, sizeArr, displ, MPI_DOUBLE, newComer->vz, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Alltoallv(fxArr, sizeArr, displ, MPI_DOUBLE, newComer->fx, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(fyArr, sizeArr, displ, MPI_DOUBLE, newComer->fy, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(fzArr, sizeArr, displ, MPI_DOUBLE, newComer->fz, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
    
    MPI_Alltoallv(massArr, sizeArr, displ, MPI_DOUBLE, newComer->mass, transSizeArr, transDispl, MPI_DOUBLE, MPI_COMM_WORLD);
   
   //each octant has its own newComer size, ncSize 
   newComer->size = newComer->used = ncSize;
}

// Function to print center of mass and average velocity
/*modified from:
  Author: Andrew Sherman, Yale University

  Date: 2/23/2013

*/
void output(int ts, Body* myOct) {
  int thisbody;
  double cmassx, cmassy, cmassz, tmass, tvelx, tvely, tvelz;
  double cmassxAll, cmassyAll, cmasszAll, tmassAll, tvelxAll, tvelyAll, tvelzAll;
  if(rank == 0){
      if (ts==0) printf("\n\nInitial Conditions (time = 0.0):\n");
      else printf("\n\nConditions after timestep %d (time = %f):\n",ts,ts*dt);
  }

  cmassx = 0.;
  cmassy = 0.;
  cmassz = 0.;
  tmass = 0.;
  tvelx = 0.;
  tvely = 0.;
  tvelz = 0.;
  
  cmassxAll = 0.;
  cmassyAll = 0.;
  cmasszAll = 0.;
  tmassAll = 0.;
  tvelxAll = 0.;
  tvelyAll = 0.;
  tvelzAll = 0.;
  for (thisbody=0; thisbody < myOct->size; thisbody++) {
    cmassx += myOct->mass[thisbody]*myOct->x[thisbody];
    cmassy += myOct->mass[thisbody]*myOct->y[thisbody];
    cmassz += myOct->mass[thisbody]*myOct->z[thisbody];
    tmass  += myOct->mass[thisbody];
    tvelx += myOct->vx[thisbody];
    tvely += myOct->vy[thisbody];
    tvelz += myOct->vz[thisbody];
  }

  barrier(); 
  //FIXME reduce problem
  MPI_Reduce(&cmassx, &cmassxAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&cmassy, &cmassyAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&cmassz, &cmasszAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&tmass, &tmassAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&tvelx, &tvelxAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&tvely, &tvelyAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
  MPI_Reduce(&tvelz, &tvelzAll, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

  if(rank == 0){
      printf("\n     Center of Mass:   (%e, %e, %e)\n", cmassxAll/tmassAll, cmassyAll/tmassAll, cmasszAll/tmassAll);
      printf(  "     Average Velocity: (%e, %e, %e)\n", tvelxAll/N, tvelyAll/N, tvelzAll/N);
  }
}

