#include "matmul.h"

//#include "timing.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
//avoid implicit declaration
void timing(double* wcTime, double* cpuTime);

void printMat(int size, double** mat){
    int i,j;

    for(i=0; i<size; i++){
        for(j=0; j<size; j++){
            printf(" %.4f ",(*mat)[i*size+j]);
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

void initAB(int sizeAB, double** A, double** B){
   int i; 
   initRowBlk(sizeAB, A);
   initColBlk(sizeAB, B);
   srand(12345); // Use a standard seed value for reproducibility

    // This assumes A is stored by rows, and B is stored by columns. Other storage schemes are permitted
    for (i=0; i<sizeAB; i++){
         (*A)[i] = ((double) rand()/(double)RAND_MAX);
    }
    for (i=0; i<sizeAB; i++){
         (*B)[i] = ((double) rand()/(double)RAND_MAX);
    }
}

void initRowBlk(int sizeA, double** A)
{
    *A = (double *) calloc(sizeA,sizeof(double));//result mat
}

void sliceBlk(int* offset, int sizeAB, double** AB, double** mat){

    memcpy(*AB, (*mat)+(*offset), sizeAB*sizeof(double));
    (*offset)=(*offset)+sizeAB;
}

void initColBlk(int sizeB, double** B)
{
    
    *B = (double *) calloc(sizeB,sizeof(double));//upper triangular mat
    
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
  int i,ii, j,jj, k;
  int iA,iB, jB, iC;
  int a1,b1;
  double wctime0, wctime1, cputime;

  timing(&wctime0, &cputime);
    
  //FIXME
  //this algorithm is modified from the serial version
  //a1 and b1 are the first element of a number series [1,2,3,4,5...n]
  //the a1 and b1 equation are drived from arithmatic sequence equation:
  //an = a1 + (n-1)*d; Sn = n*(a1+an)/2; 
  // n -> blockSize; d -> 1; Sn -> sizeAB;
  a1 = (2*sizeA/blockSize - blockSize + 1)/2;
  b1 = (2*sizeB/blockSize - blockSize + 1)/2;
  
  //iA = rankB;//col entry in each block depends on the tag of B mat
  iC = rankB*blockSize;//row entry offset based on the blocksize
  ii=0;
  for (i=a1; i<a1+blockSize; i++, ii++) {
      iA = ii*(2*a1+ii-1)/2;
      jB = 0;
      jj=0;
      for (j=b1; j<b1+blockSize; j++,jj++,jB++) {
          iB = jj*(2*b1+jj-1)/2;
          (*C)[iC+jB] = 0.;
          for (k=0; k<MIN(i,j); k++) {
              (*C)[iC+jB] += (*A)[k+iA] * (*B)[k+iB]; 
          }
      }
      iC+=cCol;//line up the offsets for next result if available
  }
  //debugPrints(C, blockSize * 3, rankB);
  timing(&wctime1, &cputime);
  return(wctime1 - wctime0);
}
