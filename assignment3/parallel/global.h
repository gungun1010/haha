#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NUM_PROCESSORS 8
#define LOCATION_TAG 0xB000
#define VELOCITY_TAG 0xC000
#define ROOT 0

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

typedef struct {
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *mass;

    size_t used;
    size_t size;
} Body;
