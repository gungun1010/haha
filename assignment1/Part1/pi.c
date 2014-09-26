/*
 *  * Author: Leon Yu, netID: lly6
 *   * department: Electrical Engineering, computer track
 *   
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SLICES 1000000000 //1 billion slices

//avoid implicit declaration
void timing(double* wcTime, double* cpuTime);

/* simple pi calculation, 
 * this function doesn't involve any unrolling 
 */
void simplePi(double* piRef, double* delta_xRef){
    int i;  //i is the index for each slice
    double x,sum;
    
    sum = 0; //init sum to 0
    //algorithem provided by cs424 assignemnt 1 doc
    for (i=0; i < SLICES; i++) {
        x = ((double)i + 0.5) * (*delta_xRef);
        sum += (1.0 / (1.0 + x * x));
    }
    *piRef = 4.0 * sum * (*delta_xRef); 
}

/*this is also a pi calculation
 * but this function unrolled for loop 10 times, as I calculated that div op takes 10 cycles/stages to complete
 * the iterations are bottlenecked by this lengthy div op
 */
void unrollPi_10(double* piRef, double* delta_xRef){
    int i;  //i is the index for each slice
    double x[10];
    double sum[10];
    double totalSum;
    
    totalSum = 0; //init sum to 0
    memset(sum,0,10*sizeof(double));
  
    //algorithem provided by cs424 assignemnt 1 doc
    for (i=0; i < SLICES; i+=10) {
        x[0] = ((double)i + 0.5) * (*delta_xRef);
        x[1] = ((double)i+1 + 0.5) * (*delta_xRef);
        x[2] = ((double)i+2 + 0.5) * (*delta_xRef);
        x[3] = ((double)i+3 + 0.5) * (*delta_xRef);
        x[4] = ((double)i+4 + 0.5) * (*delta_xRef);
        x[5] = ((double)i+5 + 0.5) * (*delta_xRef);
        x[6] = ((double)i+6 + 0.5) * (*delta_xRef);
        x[7] = ((double)i+7 + 0.5) * (*delta_xRef);
        x[8] = ((double)i+8 + 0.5) * (*delta_xRef);
        x[9] = ((double)i+9 + 0.5) * (*delta_xRef);
        sum[0] += (1.0 / (1.0 + x[0] * x[0]));
        sum[1] += (1.0 / (1.0 + x[1] * x[1]));
        sum[2] += (1.0 / (1.0 + x[2] * x[2]));
        sum[3] += (1.0 / (1.0 + x[3] * x[3]));
        sum[4] += (1.0 / (1.0 + x[4] * x[4]));
        sum[5] += (1.0 / (1.0 + x[5] * x[5]));
        sum[6] += (1.0 / (1.0 + x[6] * x[6]));
        sum[7] += (1.0 / (1.0 + x[7] * x[7]));
        sum[8] += (1.0 / (1.0 + x[8] * x[8]));
        sum[9] += (1.0 / (1.0 + x[9] * x[9]));
    }

    totalSum = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7]+sum[8]+sum[9];

   *piRef = 4.0 * totalSum * (*delta_xRef); 
}

void unrollPi_20(double* piRef, double* delta_xRef){
    int i;  //i is the index for each slice
    double x[5];
    double sum[5];
    double totalSum;
    
    totalSum = 0; //init sum to 0
    memset(sum,0,5*sizeof(double));
  
    //algorithem provided by cs424 assignemnt 1 doc
    for (i=0; i < SLICES; i+=5) {
        x[0] = ((double)i + 0.5) * (*delta_xRef);
        x[1] = ((double)i+1 + 0.5) * (*delta_xRef);
        x[2] = ((double)i+2 + 0.5) * (*delta_xRef);
        x[3] = ((double)i+3 + 0.5) * (*delta_xRef);
        x[4] = ((double)i+4 + 0.5) * (*delta_xRef);
        sum[0] += (1.0 / (1.0 + x[0] * x[0]));
        sum[1] += (1.0 / (1.0 + x[1] * x[1]));
        sum[2] += (1.0 / (1.0 + x[2] * x[2]));
        sum[3] += (1.0 / (1.0 + x[3] * x[3]));
        sum[4] += (1.0 / (1.0 + x[4] * x[4]));
    }
    totalSum = sum[0]+sum[1]+sum[2]+sum[3]+sum[4];

   *piRef = 4.0 * totalSum * (*delta_xRef); 
}


int main(){
    //declare variable for pi calculation
    double pi;
    double delta_x = 1.0/SLICES; //step size for each slice
    double wcStart,cpuStart,wcEnd,cpuEnd; //declare the time stamp variables

    //referenced from ~ahs3/cpsc424/utils/timing/timing.o 
    timing(&wcStart,&cpuStart);//time stamp start
    simplePi(&pi, &delta_x);
    timing(&wcEnd, &cpuEnd);
    //stdout
    printf("pi = %.10lf\n",pi);
    printf("wall Clock Time for simple calculation = %f\n",(wcEnd - wcStart));
  
    //referenced from ~ahs3/cpsc424/utils/timing/timing.o 
    timing(&wcStart,&cpuStart);//time stamp start
    unrollPi_20(&pi, &delta_x); 
    timing(&wcEnd,&cpuEnd);//time stamp end
    //stdout
    printf("pi = %.10lf\n",pi);
    printf("wall Clock Time for unrolling = %f\n",(wcEnd - wcStart));

    return 0;
}
