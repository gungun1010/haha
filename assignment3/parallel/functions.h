#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "/home/fas/hpcprog/ahs3/cpsc424/utils/timing/timing.h"

void printOct(Body *a,int i);

void initBody(Body *a, size_t initialSize);

void insertBody(Body *a, int i);

void freeBody(Body *a);

void output(int ts);

void force(int body1, int body2, double *deltaf);

void init(int* argcRef, char ***argvRef, int* procNumRef, int* rankRef);

void initData();

void boardcastConsts();

void scatOctants();
