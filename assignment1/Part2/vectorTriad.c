/*
 * Author: Leon Yu, netID: lly6
 * department: Electrical Engineering, computer track
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dummy.h"

#define EXPO_START 3
#define EXPO_END 24
#define BASE 2.1
#define ARRAY_VAL 1.0

int size;//global array size, so that we dont have to pass it around

//declare here to avoid implicit warning
void timing(double* wallTime, double* cpuTime);

//initialize array by using calloc(), each element is 0.0
void arrayInit(double** arrRef){
    int i;
    //allocate memory here
    *arrRef = (double*)malloc(size * sizeof(double));
    
    //init array
    for(i = 0; i < size; i++){
        (*arrRef)[i] = ARRAY_VAL;
    }
}

//benchmark function
//the outter layer is taken from assignment 1 doc
//my code starts as indicated
void vtBenchmark(double** aRef, double** bRef, double** cRef, double** dRef){
    unsigned long int repeat=1;//use unsigned long here just in case
    int r;
    double runtime=0.;
    double wcs, ct, wce;
    int i; //size indicator
    
    for(;runtime < 1; repeat*=2){//the forever loop will break out when runtime >= 1 second
        timing(&wcs, &ct);
        for(r=0; r<repeat; ++r) {
            /* BELOW ARE THE BENCHMARK LOOP */
            for(i=0; i<size; ++i){
                (*aRef)[i] = (*bRef)[i] + (*cRef)[i] * (*dRef)[i];
            }
            if ((*aRef)[size>>1] < 0.) 
                dummy(aRef); // fools the compiler
        } 
        timing(&wce, &ct);
        runtime = wce - wcs;
    }
    printf("%d      %.10lf      %d\n",size, runtime,repeat);
    repeat/=2;
}

int main(){
    double* a;
    double* b;
    double* c;
    double* d;
    int k;
   
    
    printf("N       runtime        repeat\n"); 
    for(k = EXPO_START; k <= EXPO_END; k++){ 
        size = floor(pow(BASE,k));
        
        //allocate and init arraies
        arrayInit(&a);
        arrayInit(&b);
        arrayInit(&c);
        arrayInit(&d);

        vtBenchmark(&a,&b,&c,&d);
        
        //de-allocate arries
        free(a);
        free(b);
        free(c);
        free(d);
    }
    
    
    return 0;
}
