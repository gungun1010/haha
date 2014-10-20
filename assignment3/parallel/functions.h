#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"
void barrier();

void printOct(Body *a,int i);

void initBody(Body *a, size_t initialSize);

void insertBody(Body *a, int i);

void freeBody(Body *a);

void freeOctants(Body **octRef);

void output(int ts);

void force(int body1, int body2, double *deltaf);

void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef);

void initData();

void freeInitData();

void boardcastConsts();

void sliceOctants(Body** oct);

void prepScat(Body** oct);

void scatOctants(Body* myOct);

void recvBodies(Body* myOct);

void pointToAxis(double x, double y, double z, double** duAxis);

void checkNeighbor(int a, int b, int c, double** duPlane, Body** wildCardsTo, Body* myOct, int i);

void checkOppose(int a, double duOrigin, Body** wildCardsout, Body* myOct, int i);

void freeBuffer();

void checkKeen(int a, int b, int c,double** duAxis, Body** wildCardsTo, Body* myOct, int i);

void estimateDU(Body* myOct, Body** wildCardsTo);

void prepScatWildcards(Body** wildCardsTo);

void exchangeCards(Body* myWildCards);
 
