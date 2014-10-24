#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
//MPI_Barrier() wrapper
void barrier();

//print bodies in the octant
void printOct(Body *a,int i);

//initializing Bodies in the octant
void initBody(Body *a, size_t initialSize);

//add a body to the octant
void insertBody(Body *a, int i);

//remove a body in the octant
void removeBody(Body *a, int index);

//de-allocate memory in the octant
void freeBody(Body *a);

//de-allocate memory in the my wildcards
void freeWildCards(Body *a);

//de-allocate memory in the wildcards' decks
void freeCardsDeck(Body **a);

//de-allocate memory for new comers
void freeNewcomer(Body *a);

//de-allocate everything in the space
void freeOctants(Body **octRef);

//output function taken from Dr.Sherman
void output(int ts);

//force calculation function
void force(int body1, int body2, Body* myOct, Body* myWildCards, double **deltaf);

//MPI essentials for intilaliztion
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef);

//read in data from file
void initData();

//free buffers for initial data
void freeInitData();

//MPI_Bcast all the constants
void boardcastConsts();

//split bodies into octants
void sliceOctants(Body** oct);

//prepare buffers to scat bodies to their octants
void prepScat(Body** oct);

//scat prepared buffers to their octants
void scatOctants(Body* myOct);

//calculates a DU from a point to an axis
void pointToAxis(double x, double y, double z, double** duAxis);

//check if the body is within 5DU of its neighbor octant
void checkNeighbor(int a, int b, int c, double** duPlane, Body** wildCardsTo, Body* myOct, int i);

//check if the body is within 5DU of the  octant sharing origin
void checkOppose(int a, double duOrigin, Body** wildCardsout, Body* myOct, int i);

//free up MPI buffers
void freeBuffer(int stage);


//check if the body is within 5DU of the  octant sharing an edge
void checkKeen(int a, int b, int c,double** duAxis, Body** wildCardsTo, Body* myOct, int i);

//function to estimate DU
void estimateDU(Body* myOct, Body** wildCardsTo);

void prepScatWildcards(Body** wildCardsTo);

void exchangeCards(Body* myWildCards);

void calcForce(Body* myOct, Body* myWildCards); 

int findOwner(Body* myOct, int i);

void addToOwner(Body* Owner, Body* myOct, int i);

void updateOwner(Body** oct, Body* myOct);

void prepScatNewcomer(Body** oct);

void exchangeNewcomer(Body* newComer);

void insertNewcomer(Body* a, Body* newComer, int i);

void welcomeNewcomer(Body* myOct, Body* newComer);
