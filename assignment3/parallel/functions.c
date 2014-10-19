#include "functions.h"
#include "mpi.h"

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
    //scat initOctSize of each octant for the owner
    MPI_Scatter(initOctSize, 1, MPI_INT, &myInitOctSize,1, MPI_INT, ROOT, MPI_COMM_WORLD);
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

    initOctSize = calloc(procNum, sizeof(int));
    
    //push each octant size into an array then scatter to the owner 
    for(i=0;i<procNum;i++){
        //printOct(&(*oct)[i],i);
        initOctSize[i]=(*oct)[i].size;
        //printf("%d \n",initOctSize[i]);
    }
    
}

void scatOctants(Body** oct){
   int i,j;
   MPI_Request request;

   for(i=0;i<procNum;i++){
       //distribute coordinates
        MPI_Isend((*oct)[i].x, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_X, MPI_COMM_WORLD,&request);
        MPI_Isend((*oct)[i].y, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_Y, MPI_COMM_WORLD,&request);
        MPI_Isend((*oct)[i].z, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_Z, MPI_COMM_WORLD,&request);

        //distribute velocity
        MPI_Isend((*oct)[i].vx, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_VX, MPI_COMM_WORLD,&request);
        MPI_Isend((*oct)[i].vy, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_VY, MPI_COMM_WORLD,&request);
        MPI_Isend((*oct)[i].vz, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_VZ, MPI_COMM_WORLD,&request);
        
        //distribute mass
        MPI_Isend((*oct)[i].mass, (*oct)[i].size, MPI_DOUBLE, i, i+TAG_MASS, MPI_COMM_WORLD,&request);
   }
    
}

void recvBodies(Body* myOct){
   //we can hardcode 7 here since the parameters are fixed to [x,y,z], [vx,vy,vz], mass
   MPI_Request request[7];
   MPI_Status status[7];
   
   //allocate coordinates for bodies in the octant 
   myOct->x = (double *) calloc(myInitOctSize, sizeof(double)); 
   myOct->y = (double *) calloc(myInitOctSize, sizeof(double)); 
   myOct->z = (double *) calloc(myInitOctSize, sizeof(double)); 
   
   //allocate velocity for bodies in the octant
   myOct->vx = (double *) calloc(myInitOctSize, sizeof(double)); 
   myOct->vy = (double *) calloc(myInitOctSize, sizeof(double)); 
   myOct->vz = (double *) calloc(myInitOctSize, sizeof(double)); 
    
   //allocate masses for bodies in the octant
   myOct->mass = (double *) calloc(myInitOctSize, sizeof(double)); 

   //wait to get coordinates
   MPI_Irecv(myOct->x, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_X, MPI_COMM_WORLD,&request[0]);
   MPI_Irecv(myOct->y, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_Y, MPI_COMM_WORLD,&request[1]);
   MPI_Irecv(myOct->z, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_Z, MPI_COMM_WORLD,&request[2]);
   
   //wait to get velocity
   MPI_Irecv(myOct->vx, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_VX, MPI_COMM_WORLD,&request[3]);
   MPI_Irecv(myOct->vy, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_VY, MPI_COMM_WORLD,&request[4]);
   MPI_Irecv(myOct->vz, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_VZ, MPI_COMM_WORLD,&request[5]);

   //wait to get mass
   MPI_Irecv(myOct->mass, myInitOctSize, MPI_DOUBLE, ROOT, rank+TAG_MASS, MPI_COMM_WORLD,&request[6]);
   
   //wait all to finish
   MPI_Waitall(7,request,status);
   myOct->used = myOct->size = myInitOctSize;
}

void pointToAxis(double x, double y, double z, double** duAxis){
    double dx,dy,dz;

    dx = x-x;
    dy = y-0.;
    dz = z-0.;
    (*duAxis)[0]=sqrt(dx*dx + dy*dy + dz*dz);
    
    dx = x-0.;
    dy = y-y;
    dz = z-0.;
    (*duAxis)[1]=sqrt(dx*dx + dy*dy + dz*dz);
    
    dx = x-0.;
    dy = y-0.;
    dz = z-z;
    (*duAxis)[2]=sqrt(dx*dx + dy*dy + dz*dz);
}

void estimateDU(Body* myOct, Body* wildCard){
    int i;
    double x,y,z;
    double duXY;
    double duYZ;
    double duZX;
    double* duAxis;
    double duOrigin;
    
    //still dynamically allocate even though i know it is always gonna be 3 axises 
    //duAxis[0] -----> x axis DU
    //duAxis[1] -----> y axis DU
    //duAxis[2] -----> z axis DU
    duAxis = calloc(3,sizeof(double));
    
    for(i=0;i<myOct->size;i++){
        x = myOct->x[i];
        y = myOct->y[i];
        z = myOct->z[i];
        
        //distance between the point to planes
        //only send location and mass
        duXY = fabs(z);
        duYZ = fabs(x);
        duZX = fabs(y);   
        
        //distance between the point to axis
        pointToAxis(x,y,z,&duAxis);

        //distance between the point to origin
        duOrigin = sqrt(x*x + y*y + z*z);

        switch (rank){
            case 0://neighbor 1, 3, 4, oppose 6
                if(duXY <= 5.) sendWildCard(4);
                if(duYZ <= 5.) sendWildCard(3);
                if(duZX <= 5.) sendWildCard(1);
            case 1://neighbor 0, 2, 5, oppose 7
            case 2://neighbor 1, 3, 6, oppose 4 
            case 3://neighbor 0, 7, 2, oppose 5
            case 4://neighbor 0, 5, 7, oppose 2
            case 5://neighbor 1, 4, 6, oppose 3
            case 6://neighbor 2, 5, 7, oppose 0
            case 7://neighbor 3, 4, 6, oppose 1
    }
        
}







