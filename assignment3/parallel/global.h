#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NUM_PROCESSORS 8

#define ROOT 0
#define DU_THRES 5.
#define INIT_STAGE 0xA
#define WILDCARD_STAGE 0xB
#define NEWCOMER_STAGE 0xC

//data type for each body, contains its location, velocity, and mass
typedef struct {
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *mass;
    
    //these are the force calc results, they never scats
    double *fx;     // X force array for N bodies
    double *fy;     // Y velocity array for N bodies
    double *fz;     // Z velocity array for N bodies

    //these two are for initial distribution
    size_t used;
    size_t size;
} Body;

// Globals so that we don't need to pass them to functions
double dt;      //delta t for [t_k, t_k+1]
double dtBy2;   //delta t divided by 2
int N;          // N bodies
int K;          // K time steps
double *mass;   // mass array storing N bodies' mass
double *x;      // x array storing x axis coordinates for N bodies
double *y;      // y array storing y axis coordinates for N bodies
double *z;      // z array storing z axis coordinates for N bodies
double *vx;     // X velocity array for N bodies
double *vy;     // Y velocity array for N bodies
double *vz;     // Z velocity array for N bodies
double *fx;     // X force array for N bodies
double *fy;     // Y velocity array for N bodies
double *fz;     // Z velocity array for N bodies
int procNum;    // process number in the communicator
int rank;       // rank of proc in the communicator
int *sizeArr, *displ; //array of initial size of each octant, and displacement for the array
int *transSizeArr, *transDispl; // buffer to store transposed size and displacement for alltoall operation
int mySize; //size of each octant, used in each process   

//buffers for collective comm.
double *xArr, *yArr, *zArr, *vxArr, *vyArr, *vzArr, *fxArr, *fyArr, *fzArr, *massArr;
