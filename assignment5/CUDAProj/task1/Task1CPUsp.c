#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cpu_matrixmult(float *a,float *b, float **c, int n, int m, int p) {

  int k,i,j;
  float r;

  /*kij*/
  for (k=0; k<p; k++){
      for (i=0; i<n; i++){
          r = a[i*p + k];
          for(j=0; j<m; j++){
              (*c)[i*m+j] += r * b[k*m+j];
          }
      }
  }
}


int main(int argc, char *argv[]) {

  int i, j; // loop counters

  int gpucount = 0; // Count of available GPUs
  int Grid_Dim = 1; //Grid dimension, x and y, square
  int Block_Dim = 1; //Block dimension, x and y, square

  int n,m,p; // matrix dimension
  float *a,*b,*c;
  int size_a, size_b, size_c; // number of bytes in arrays

  float elapsed_time_ms; // which is applicable for asynchronous code also

  // --------------------SET PARAMETERS AND DATA -----------------------


  if (argc<4) {
    printf("Usage: Task1GPUsp <n> <m> <p>\n");
    exit (-1);
  }

  n = atoi(argv[1]);
  m = atoi(argv[2]);
  p = atoi(argv[3]);
  
  
  printf("A Matrix Dimension = %dx%d\n",n,p);
  printf("B Matrix Dimension = %dx%d\n",p,m);
  printf("C Matrix Dimension = %dx%d\n",n,m);

  size_a = n * p * sizeof(float); // number of bytes in total in arrays
  size_b = p * m * sizeof(float); // number of bytes in total in arrays
  size_c = n * m * sizeof(float); // number of bytes in total in arrays

  a = (float*) malloc(size_a); // dynamically allocated memory for arrays on host
  b = (float*) malloc(size_b);
  c = (float*) malloc(size_c); // results from CPU

  srand(12345);
  //int p = n; //Used here only to illustrate proper initialization for non-square case
  printf ("a\n");
  for(i=0;i < n;i++){
    for(j=0;j < p;j++) {
      a[i * p + j] = (float) rand() / (float) RAND_MAX;
      //a[i * p + j] = (float) (i+j);
      printf("%.2f  ", a[i * p + j]);
    }
    printf("\n");
  }

  printf("b\n");
  for(i=0;i < p;i++){
    for(j=0;j < m;j++) {
      b[i * m + j] = (float) rand() / (float) RAND_MAX;
      //b[i * m + j] = (float) (i+j);
      printf("%.2f  ", b[i * m + j]);
    }
    printf("\n");
  }

  // ------------- COMPUTATION DONE ON HOST CPU ----------------------------

  cpu_matrixmult(a,b, &c, n, m, p); //kij implementation)

  printf("c\n");
  for(i=0;i < n;i++){
    for(j=0;j < m;j++) {
      printf("%.2f  ", c[i * m + j]);
    }
    printf("\n");
  }

  free(a);
  free(b);
  free(c);
  return 0;
}
