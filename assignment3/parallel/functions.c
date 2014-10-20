#include "functions.h"
#include "mpi.h"

void barrier(){
    MPI_Barrier(MPI_COMM_WORLD);
}

void printOct(Body *a, int i){
    int j;
    printf("octant %d *********************************************************\n",i);
    
    for(j=0;j<a->used;j++){
        printf("coordinate: [%.3le, %.3le, %.3le] ",a->x[j],a->y[j],a->z[j]);
        printf("velocity: [%.3le, %.3le, %.3le] ",a->vx[j],a->vy[j],a->vz[j]);
        printf("mass: %.3le\n",a->mass[j]);
    }
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

void freeBody(Body *a) {
      free(a->x);
      free(a->y);
      free(a->z);
      free(a->vx);
      free(a->vy);
      free(a->vz);
      free(a->mass);

      a->used = a->size = 0;
}

void freeOctants(Body** octRef){
    int i;
    for(i=0;i<procNum;i++){
        freeBody(&(*octRef)[i]);
    }
}

void freeBuffer(){
    free(massArr);
    free(xArr);
    free(yArr);
    free(zArr);
    free(vxArr);
    free(vyArr);
    free(vzArr);
    free(displ);
    free(sizeArr);
}
// Function to print center of mass and average velocity
/*
  Author: Andrew Sherman, Yale University

  Date: 2/23/2013

*/
void output(int ts) {
  int thisbody;
  double cmassx, cmassy, cmassz, tmass, tvelx, tvely, tvelz;

  if (ts==0) printf("\n\nInitial Conditions (time = 0.0):\n");
  else printf("\n\nConditions after timestep %d (time = %f):\n",ts,ts*dt);

  cmassx = 0.;
  cmassy = 0.;
  cmassz = 0.;
  tmass = 0.;
  tvelx = 0.;
  tvely = 0.;
  tvelz = 0.;
  for (thisbody=0; thisbody<N; thisbody++) {
    cmassx += mass[thisbody]*x[thisbody];
    cmassy += mass[thisbody]*y[thisbody];
    cmassz += mass[thisbody]*z[thisbody];
    tmass  += mass[thisbody];
    tvelx += vx[thisbody];
    tvely += vy[thisbody];
    tvelz += vz[thisbody];
  }
  printf("\n     Center of Mass:   (%e, %e, %e)\n", cmassx/tmass, cmassy/tmass, cmassz/tmass);
  printf(  "     Average Velocity: (%e, %e, %e)\n", tvelx/N, tvely/N, tvelz/N);
}

// Function to compute the forces between a pair of bodies
/*
  Author: Andrew Sherman, Yale University

  Date: 2/23/2013

*/
void force(int body1, int body2, double *deltaf) {
  double gmmr3, r, r2, dx, dy, dz;
  double G=1.0;
  
  dx = x[body2] - x[body1];
  dy = y[body2] - y[body1];
  dz = z[body2] - z[body1];
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  if (r<=5.) {
    gmmr3 = G*mass[body1]*mass[body2]/(r2*r);
    deltaf[0] = gmmr3 * dx;
    deltaf[1] = gmmr3 * dy;
    deltaf[2] = gmmr3 * dz;
  }
  else{
    deltaf[0] = 0.;
    deltaf[1] = 0.;
    deltaf[2] = 0.;
  }
}

//MPI initialization function
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef){
    
    MPI_Init(argcRef, argvRef);

    MPI_Comm_size(MPI_COMM_WORLD, procNumRef); // Get # of processes from MPI commnad
    MPI_Comm_rank(MPI_COMM_WORLD, rankRef); // Which process am I?
}

//function to init all data, using master
/*
  Author: Andrew Sherman, Yale University

  Date: 2/23/2013

*/
void initData(){
    int i;
    
    // Read basic inputs   
    scanf("%d\n",&N); //first line of file
    scanf("%d\n",&K); //2nd line of file
    scanf("%le\n",&dt); //thrid line of file
     
    dtBy2 = dt/2.;
    
    // Allocate arrays
    mass = calloc(N, sizeof(double)); // mass[i] is mass of body i
    x = calloc(N, sizeof(double)); // x[i] is x position of body i
    y = calloc(N, sizeof(double)); // y[i] is y position of body i
    z = calloc(N, sizeof(double)); // y[i] is z position of body i
    vx = calloc(N, sizeof(double)); // vx[i] is x velocity of body i
    vy = calloc(N, sizeof(double)); // vy[i] is y velocity of body i
    vz = calloc(N, sizeof(double)); // vz[i] is z velocity of body i
      
    // Read initial conditions
    for (i=0; i<N; i++) scanf("%le\n",&mass[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&x[i],&y[i],&z[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&vx[i],&vy[i],&vz[i]);
}

void freeInitData(){
    free(mass);
    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
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

void sliceOctants(Body** oct){
    int i;
    
    *oct = (Body *) calloc(procNum, sizeof(Body));

    for(i=0;i<procNum;i++){
        initBody(&(*oct)[i], 1);  // initially 1 elements     
    }
   
    for(i=0;i<N; i++){
        if(x[i]>=0. && y[i]>=0. && z[i]>=0.)//0
            insertBody(&(*oct)[0],i);
        if(x[i]>=0. && y[i]<0. && z[i]>=0.)//1
            insertBody(&(*oct)[1],i);
        if(x[i]<0. && y[i]<0. && z[i]>=0.)//2
            insertBody(&(*oct)[2],i);
        if(x[i]<0. && y[i]>=0. && z[i]>=0.)//3
            insertBody(&(*oct)[3],i);
        if(x[i]>=0. && y[i]>=0. && z[i]<0.)//4
            insertBody(&(*oct)[4],i);
        if(x[i]>=0. && y[i]<0. && z[i]<0.)//5
            insertBody(&(*oct)[5],i);
        if(x[i]<0. && y[i]<0. && z[i]>=0.)//6
            insertBody(&(*oct)[6],i);
        if(x[i]<0. && y[i]>=0. && z[i]<0.)//7
            insertBody(&(*oct)[7],i);
    }

    sizeArr = calloc(procNum, sizeof(int));
    
    //push each octant size into an array then boardcast the array 
    //i do boardcast instead of scatter since the array is needed for MPI_Scatterv() 
    for(i=0;i<procNum;i++){
        sizeArr[i]=(*oct)[i].size;
    }
    
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

//FIXME unused function
void recvBodies(Body* myOct){
   //we can hardcode 7 here since the parameters are fixed to [x,y,z], [vx,vy,vz], mass
   MPI_Request request[7];
   MPI_Status status[7];
   
   //allocate coordinates for bodies in the octant 
   myOct->x = (double *) calloc(mySize, sizeof(double)); 
   myOct->y = (double *) calloc(mySize, sizeof(double)); 
   myOct->z = (double *) calloc(mySize, sizeof(double)); 
   
   //allocate velocity for bodies in the octant
   myOct->vx = (double *) calloc(mySize, sizeof(double)); 
   myOct->vy = (double *) calloc(mySize, sizeof(double)); 
   myOct->vz = (double *) calloc(mySize, sizeof(double)); 
    
   //allocate masses for bodies in the octant
   myOct->mass = (double *) calloc(mySize, sizeof(double)); 

   //wait to get coordinates
   MPI_Irecv(myOct->x, mySize, MPI_DOUBLE, ROOT, rank+TAG_X, MPI_COMM_WORLD,&request[0]);
   MPI_Irecv(myOct->y, mySize, MPI_DOUBLE, ROOT, rank+TAG_Y, MPI_COMM_WORLD,&request[1]);
   MPI_Irecv(myOct->z, mySize, MPI_DOUBLE, ROOT, rank+TAG_Z, MPI_COMM_WORLD,&request[2]);
   
   //wait to get velocity
   MPI_Irecv(myOct->vx, mySize, MPI_DOUBLE, ROOT, rank+TAG_VX, MPI_COMM_WORLD,&request[3]);
   MPI_Irecv(myOct->vy, mySize, MPI_DOUBLE, ROOT, rank+TAG_VY, MPI_COMM_WORLD,&request[4]);
   MPI_Irecv(myOct->vz, mySize, MPI_DOUBLE, ROOT, rank+TAG_VZ, MPI_COMM_WORLD,&request[5]);

   //wait to get mass
   MPI_Irecv(myOct->mass, mySize, MPI_DOUBLE, ROOT, rank+TAG_MASS, MPI_COMM_WORLD,&request[6]);
   
   //wait all to finish
   MPI_Waitall(7,request,status);
   myOct->used = myOct->size = mySize;
}

void pointToAxis(double x, double y, double z, double** duAxis){
    double dx,dy,dz;

    dx = x-x;
    dy = y-0.;
    dz = z-0.;
    (*duAxis)[0]=sqrt(dx*dx + dy*dy + dz*dz); //distance to axis x
    
    dx = x-0.;
    dy = y-y;
    dz = z-0.;
    (*duAxis)[1]=sqrt(dx*dx + dy*dy + dz*dz); //distance to axis y
    
    dx = x-0.;
    dy = y-0.;
    dz = z-z;
    (*duAxis)[2]=sqrt(dx*dx + dy*dy + dz*dz); //distance to axis z
}

void insertWildCard(Body* a, Body* myOct, int i){
      if (a->used == a->size) {
          a->size++;
          
          a->x = (double *)realloc(a->x, a->size * sizeof(double));
          a->y = (double *)realloc(a->y, a->size * sizeof(double));
          a->z = (double *)realloc(a->z, a->size * sizeof(double));
          a->mass = (double *)realloc(a->mass, a->size * sizeof(double));
      }
      
      a->x[a->used] = myOct->x[i];
      a->y[a->used] = myOct->y[i];
      a->z[a->used] = myOct->z[i];
      a->mass[a->used] = myOct->mass[i];
      a->used++;
}

void checkNeighbor(int a, int b, int c, double** duPlane, Body** wildCardsTo, Body* myOct, int i){
   
   if((*duPlane)[0] < DU_THRES){//DU to plane XY
       //add a wild card to octant a
       insertWildCard( &(*wildCardsTo)[a], myOct, i); 
   }
   if((*duPlane)[1] < DU_THRES){//DU to plane YZ
       //add a wild card to octant b 
       insertWildCard( &(*wildCardsTo)[b], myOct, i); 
   }
   if((*duPlane)[2] < DU_THRES){//DU to plane ZX
       //add a wild card to octant c
       insertWildCard( &(*wildCardsTo)[c], myOct, i);
   }
}

void checkOppose(int a, double duOrigin, Body** wildCardsTo, Body* myOct, int i){

    if(duOrigin < DU_THRES){//du to origin
       insertWildCard( &(*wildCardsTo)[a], myOct, i);
    }
}

void checkKeen(int a, int b, int c,double** duAxis, Body** wildCardsTo, Body* myOct, int i){
    if((*duAxis)[0] < DU_THRES){//DU to x axis
       //add a wild card to octant a
       insertWildCard( &(*wildCardsTo)[a], myOct, i); 
    }
    if((*duAxis)[1] < DU_THRES){//DU to y axis
       //add a wild card to octant b
       insertWildCard( &(*wildCardsTo)[b], myOct, i); 
    }
    if((*duAxis)[2] < DU_THRES){//DU to z axis
       //add a wild card to octant c
       insertWildCard( &(*wildCardsTo)[c], myOct, i); 
    }
}

void prepScatWildcards(Body** wildCardsTo){
   int i;
   int bodyNum=0;
   int offset=0;
   int size;

    //get total number of bodies in the problem
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

void estimateDU(Body* myOct, Body** wildCardsTo){
    int i;
    double x,y,z;
    double* duPlane;
    double* duAxis;
    double duOrigin;
    
    *wildCardsTo = (Body *) calloc(procNum, sizeof(Body));
    //still dynamically allocate even though i know it is always gonna be 3 axises 

    for(i=0;i<procNum;i++){
        initBody(&(*wildCardsTo)[i], 1);  // initially 1 elements     
    }

    //duPlane[0] ----> plane XY DU
    //duPlane[1] ----> plane YZ DU
    //duPlane[2] ----> plane ZX DU
    duPlane = calloc(3,sizeof(double));
    //duAxis[0] -----> x axis DU
    //duAxis[1] -----> y axis DU
    //duAxis[2] -----> z axis DU
    duAxis = calloc(3,sizeof(double));
    
    for(i=0;i<myOct->size;i++){
        x = myOct->x[i];
        y = myOct->y[i];
        z = myOct->z[i];
        
        //distance between the point to planes
        //only send coordinates and mass
        //duPlane[0] ----> plane XY DU
        //duPlane[1] ----> plane YZ DU
        //duPlane[2] ----> plane ZX DU
        duPlane[0] = fabs(z);
        duPlane[1] = fabs(x);
        duPlane[2] = fabs(y);   
        
        //distance between the point to axis
        pointToAxis(x,y,z,&duAxis);

        //distance between the point to origin
        duOrigin = sqrt(x*x + y*y + z*z);

        switch (rank){
        case 0://neighbor 4(XY), 3(YZ), 1(ZX), oppose 6, keen (sharing an edge) 5, 7, 2
            checkNeighbor(4,3,1, &duPlane, wildCardsTo, myOct, i);
            checkOppose(6, duOrigin, wildCardsTo, myOct, i);
            checkKeen(5,7,2, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 1://neighbor 5, 2, 0, oppose 7, keen 4, 6, 3
            checkNeighbor(5,2,0, &duPlane, wildCardsTo, myOct, i);
            checkOppose(7, duOrigin, wildCardsTo, myOct, i);
            checkKeen(4,6,3, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 2://neighbor 6, 1, 3, oppose 4, keen 7, 5, 0 
            checkNeighbor(6,1,3, &duPlane, wildCardsTo, myOct, i);
            checkOppose(4, duOrigin, wildCardsTo, myOct, i);
            checkKeen(7,5,0, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 3://neighbor 7, 0, 2, oppose 5, keen 6, 4, 1
            checkNeighbor(7,0,2, &duPlane, wildCardsTo, myOct, i);
            checkOppose(5, duOrigin, wildCardsTo, myOct, i);
            checkKeen(6,4,1, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 4://neighbor 0, 7, 5,  oppose 2, keen 1, 3, 6
            checkNeighbor(0,7,5, &duPlane, wildCardsTo, myOct, i);
            checkOppose(2, duOrigin, wildCardsTo, myOct, i);
            checkKeen(1,3,6, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 5://neighbor 1, 6, 4,  oppose 3, keen 0, 2, 7
            checkNeighbor(1,6,4, &duPlane, wildCardsTo, myOct, i);
            checkOppose(3, duOrigin, wildCardsTo, myOct, i);
            checkKeen(0,2,7, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 6://neighbor 2, 5, 7, oppose 0, keen 3, 1, 4
            checkNeighbor(2,5,7, &duPlane, wildCardsTo, myOct, i);
            checkOppose(0, duOrigin, wildCardsTo, myOct, i);
            checkKeen(3,1,4, &duAxis, wildCardsTo, myOct, i); 
            break;
        case 7://neighbor 3, 4, 6, oppose 1, keen 2, 0, 5
            checkNeighbor(3,4,6, &duPlane, wildCardsTo, myOct, i);
            checkOppose(1, duOrigin, wildCardsTo, myOct, i);
            checkKeen(2,0,5, &duAxis, wildCardsTo, myOct, i); 
            break;
        default:
            break;
        }
    }
}

