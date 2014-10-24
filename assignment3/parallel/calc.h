#include "comm.h"

//print bodies in the octant
void printOct(Body *a,int i);

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


//force calculation function
void force(int body1, int body2, Body* myOct, Body* myWildCards, double **deltaf);

//read in data from file
void initData();

//free buffers for initial data
void freeInitData();

//split bodies into octants
void sliceOctants(Body** oct);

//calculates a DU from a point to an axis
void pointToAxis(double x, double y, double z, double** duAxis);

void insertWildCard(Body* a, Body* myOct, int i);
//check if the body is within 5DU of its neighbor octant
void checkNeighbor(int a, int b, int c, double** duPlane, Body** wildCardsTo, Body* myOct, int i);

//check if the body is within 5DU of the  octant sharing origin
void checkOppose(int a, double duOrigin, Body** wildCardsout, Body* myOct, int i);

//check if the body is within 5DU of the  octant sharing an edge
void checkKeen(int a, int b, int c,double** duAxis, Body** wildCardsTo, Body* myOct, int i);

//function to estimate DU
void estimateDU(Body* myOct, Body** wildCardsTo);

//calculate force
void calcForce(Body* myOct, Body* myWildCards); 

//find who is my new owner
int findOwner(Body* myOct, int i);

//add new body to the owner list
void addToOwner(Body* Owner, Body* myOct, int i);

//update the owner
void updateOwner(Body** oct, Body* myOct);

//insert new comers into my Octant
void insertNewcomer(Body* a, Body* newComer, int i);

void welcomeNewcomer(Body* myOct, Body* newComer);
