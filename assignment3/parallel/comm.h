#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"

//MPI_Barrier() wrapper
void barrier();

//initializing Bodies in the octant
void initBody(Body *a, size_t initialSize);

//add a body to the octant
void insertBody(Body *a, int i);

//remove a body in the octant
void removeBody(Body *a, int index);


//MPI essentials for intilaliztion
void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef);

//MPI_Bcast all the constants
void boardcastConsts();

//prepare buffers to scat bodies to their octants
void prepScat(Body** oct);

//scat prepared buffers to their octants
void scatOctants(Body* myOct);

//free up MPI buffers
void freeBuffer(int stage);

//prepare buffer to scat wild cards
void prepScatWildcards(Body** wildCardsTo);

//all to all exchange cards
void exchangeCards(Body* myWildCards);

//prepare scat new comer list
void prepScatNewcomer(Body** oct);

//all to all exchange new comers
void exchangeNewcomer(Body* newComer);

//output function taken from Dr.Sherman
void output(int ts, Body* myOct);
