#include "matmul.h"

//#include "timing.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
//avoid implicit declaration
void timing(double* wcTime, double* cpuTime);

void printMat(int size, double** mat){
    int i,j;

    for(i=0; i<size; i++){
        for(j=0; j<size; j++){
            printf(" %.2f ",(*mat)[i*size+j]);
        }
        printf("\n");
    }
}

int calcSize(int rank, int blockSize)
{
    int n;
    int sizeAB;

    n = (rank+1)*blockSize;
    sizeAB = n*(1+n)/2 - (n-blockSize)*(1+n-blockSize)/2;
    return sizeAB;
}

void initRowBlk(int sizeA, int sizeC, double** A, double** C)
{

    int i;
    
    *A = (double *) calloc(sizeA,sizeof(double));//lower triangular mat
    *C = (double *) calloc(sizeC,sizeof(double));//result mat

    srand(12345);
    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeA; i++) (*A)[i] = ((double) rand()/(double)RAND_MAX);
    //for (i=0; i<sizeA; i++) (*A)[i] = 1.0;

}

void initColBlk(int sizeB, double** B)
{

    int i;
    
    *B = (double *) calloc(sizeB,sizeof(double));//upper triangular mat
    srand(12345);

    // This assumes A is stored by rows, and B is stored by columns
    for (i=0; i<sizeB; i++) (*B)[i] = ((double) rand()/(double)RAND_MAX);
    //for (i=0; i<sizeB; i++) (*B)[i] = 1.0;
    
}

void updateIndx(int* i, int procNum){
    if(*i != 0){
        (*i)--;
    }else{
        *i = procNum-1;
    }
}
void matFree(double* A, double* B, double* C)
{
    free(A);
    free(B);
    free(C);
}

void initAnC(int rank, int blockSize, int N, int* sizeAref,int* sizeCref, double** A, double** C){
    //init A row block and C row block for this process (rank 0)
    //C row block has constant size, calculated insize initRowBlk()
    *sizeAref = calcSize(rank, blockSize);
    *sizeCref = blockSize * N; 
    
    //init row block A and C based on sizeA and sizeC
    initRowBlk(*sizeAref, *sizeCref, A, C);
}

void debugPrints(double** matRef, int size, int rank){
    int i;
    printf("%d : ",rank);
    for(i = 0; i<size; i++){
        printf("%.2f ", (*matRef)[i]);
    }
    printf("\n");
}
double matmul(int rankB, int cCol, int blockSize, int sizeA, int sizeB, double** A, double** B, double** C) 
{
  int i, j, k;
  int iA, jB, iC;
  int a1,b1;
  double wctime0, wctime1, cputime;

  timing(&wctime0, &cputime);

  //this algorithm is modified from the serial version
  //a1 and b1 are the first element of a number series [1,2,3,4,5...n]
  //the a1 and b1 equation are drived from arithmatic sequence equation:
  //an = a1 + (n-1)*d; Sn = n*(a1+an)/2; 
  // n -> blockSize; d -> 1; Sn -> sizeAB;
  a1 = (2*sizeA/blockSize - blockSize + 1)/2;
  b1 = (2*sizeB/blockSize - blockSize + 1)/2;
  
  iA = rankB;//col entry in each block depends on the tag of B mat
  iC = iA*blockSize;//row entry offset based on the blocksize
  for (i=a1; i<a1+blockSize; i++, iA++) {
      jB =0;
      for (j=b1; j<b1+blockSize; j++,jB++) {
          (*C)[iC+jB] = 0.;
          for (k=0; k<MIN(i,j); k++) (*C)[iC+jB] += (*A)[k] * (*B)[k]; 
          //printf("%.1f %.1f %.1f\n",(*C)[iC+jB],(*A)[k],(*B)[k]);
      }
      iC+=cCol;//line up the offsets for next result if available
  }
  //debugPrints(C, blockSize * 3, rankB);
  timing(&wctime1, &cputime);
  return(wctime1 - wctime0);
}
