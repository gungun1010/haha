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
    fx = calloc(N, sizeof(double)); // fx[i] is x force on body i
    fy = calloc(N, sizeof(double)); // fy[i] is y force on body i
    fz = calloc(N, sizeof(double)); // fz[i] is z force on body i
      
    // Read initial conditions
    for (i=0; i<N; i++) scanf("%le\n",&mass[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&x[i],&y[i],&z[i]);
    for (i=0; i<N; i++) scanf("%le %le %le\n",&vx[i],&vy[i],&vz[i]);
    
    printf("first mass = %le\n",mass[0]);
}

void boardcastConsts(){
        
    //boardcast N, K, Dt
    MPI_Bcast(&N, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&K, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

}

void scatOctants(){
    int i;
    Body oct0;
    Body oct1;
    Body oct2;
    Body oct3;
    Body oct4;
    Body oct5;
    Body oct6;
    Body oct7;
    
    initBody(&oct0, 1);  // initially 1 elements
    initBody(&oct1, 1);  // initially 1 elements
    initBody(&oct2, 1);  // initially 1 elements
    initBody(&oct3, 1);  // initially 1 elements
    initBody(&oct4, 1);  // initially 1 elements
    initBody(&oct5, 1);  // initially 1 elements
    initBody(&oct6, 1);  // initially 1 elements
    initBody(&oct7, 1);  // initially 1 elements

    for(i=0;i<N; i++){
        if(x[i]>=0. && y[i]>=0. && z[i]>=0.)//0
            insertBody(&oct0,i);
        if(x[i]>=0. && y[i]<0. && z[i]>=0.)//1
            insertBody(&oct1,i);
        if(x[i]<0. && y[i]<0. && z[i]>=0.)//2
            insertBody(&oct2,i);
        if(x[i]<0. && y[i]>=0. && z[i]>=0.)//3
            insertBody(&oct3,i);
        if(x[i]>=0. && y[i]>=0. && z[i]<0.)//4
            insertBody(&oct4,i);
        if(x[i]>=0. && y[i]<0. && z[i]<0.)//5
            insertBody(&oct5,i);
        if(x[i]<0. && y[i]<0. && z[i]>=0.)//6
            insertBody(&oct6,i);
        if(x[i]<0. && y[i]>=0. && z[i]<0.)//7
            insertBody(&oct7,i);
    }
        

    printOct(&oct0,0);
    printOct(&oct1,1);
    printOct(&oct2,2);
    printOct(&oct3,3);
    printOct(&oct4,4);
    printOct(&oct5,5);
    printOct(&oct6,6);
    printOct(&oct7,7);
}














