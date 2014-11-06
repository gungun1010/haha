#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

void put_queue(int **vqRef, int vertex, int *qSize);

int get_queue(int **vqRef, int *qSize);

void mooresLaw(int N, int source, GNodePtr **adj_listheadRef);
