#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <math.h>

__global__ void gpu_matrixmult(double *a, double *b, double *c, int n, int m, int p) {

  int col = threadIdx.x + blockDim.x * blockIdx.x;
  int row = threadIdx.y + blockDim.y * blockIdx.y;

  int indexb = col;
  int index = row * m + col;

  if(col < m && row < n) {
    c[index] = 0.;
    for (int indexa = row*p; indexa < (row*p + p); indexa++, indexb+=m){ 
      c[index] += a[indexa]*b[indexb];
    }
  }
}


void cpu_matrixmult(double *a,double *b, double *c, int n, int m, int p) {

  int index, indexa, indexb;
  double cvalue;
  for(int col=0;col < m; col++)
    for(int row=0;row < n; row++) {
      indexb = col;
      index = row * m + col;
      cvalue = 0.;
      for (indexa = row*p; indexa < (row*p + p); indexa++, indexb+=m) 
	cvalue += a[indexa]*b[indexb];
      c[index] -= cvalue; //NOTE: This calculates the diff between CPU and GPU computations.
    }
}


int main(int argc, char *argv[]) {

  int i, j; // loop counters

  int gpucount = 0; // Count of available GPUs
  int Grid_Dim_x = 1; //Grid dimension, x
  int Grid_Dim_y = 1; //Grid dimension, y
  int Block_Dim_x = 1; //Block dimension, x
  int Block_Dim_y = 1; //Block dimension, y

  int n,m,p; // matrix dimension
  double *a,*b,*c;
  double *dev_a, *dev_b, *dev_c;
  int size_a, size_b, size_c; // number of bytes in arrays

  cudaEvent_t start, stop; // using cuda events to measure time
  float elapsed_time_ms; // which is applicable for asynchronous code also
  cudaError_t errorcode;

  // --------------------SET PARAMETERS AND DATA -----------------------

  errorcode = cudaGetDeviceCount(&gpucount);
  if (errorcode == cudaErrorNoDevice) {
    printf("No GPUs are visible\n");
    exit(-1);
  }
  //else printf("Device count = %d\n",gpucount);
  if (sizeof(argv)<8) {
    printf("Usage: Task1GPUsp <n> <m> <p> <block dim x> <block dim y> <grid dim x> <grid dim y>\n");
    exit (-1);
  }

  n = atoi(argv[1]);
  m = atoi(argv[2]);
  p = atoi(argv[3]);
  

  Block_Dim_x = atoi(argv[4]); // non-Square block, x dimension size (# of cols)
  Block_Dim_y = atoi(argv[5]); // non-Square block, y dimension size (# of rows)
  if (Block_Dim_x * Block_Dim_y > 1024) {
    printf("Error, too many threads in block\n");
    exit (-1);
  }

  Grid_Dim_x = atoi(argv[6]); // non-Square grid, x diemnsion size (# of cols)
  Grid_Dim_y = atoi(argv[7]); // non-Square grid, y dimension size (# of rows)
  if (Grid_Dim_x*Block_Dim_x < m ) {
    printf("Error, number of threads in x dimensions less than number of array elements\n");
    exit (-1);
  }

  if (Grid_Dim_y*Block_Dim_y < n) {
    printf("Error, number of threads in y dimensions less than number of array elements\n");
    exit (-1);
  }
  
  //printf("A Matrix Dimension = %dx%d\n",n,p);
  //printf("B Matrix Dimension = %dx%d\n",p,m);
  //printf("C Matrix Dimension = %dx%d\n",n,m);
  //printf("Block_x = %d Block_y = %d, Grid_x = %d Grid_y = %d\n",Block_Dim_x, Block_Dim_y,Grid_Dim_x, Grid_Dim_y);

  dim3 Grid(Grid_Dim_x, Grid_Dim_y); //Grid structure
  dim3 Block(Block_Dim_x, Block_Dim_y); //Block structure

  size_a = n * p * sizeof(double); // number of bytes in total in arrays
  size_b = p * m * sizeof(double); // number of bytes in total in arrays
  size_c = n * m * sizeof(double); // number of bytes in total in arrays

  a = (double*) malloc(size_a); // dynamically allocated memory for arrays on host
  b = (double*) malloc(size_b);
  c = (double*) malloc(size_c); // results from GPU

  srand(12345);
  //int p = n; //Used here only to illustrate proper initialization for non-square case
  
  //printf ("a\n");
  for(i=0;i < n;i++){
    for(j=0;j < p;j++) {
      a[i * p + j] = (double) rand() / (double) RAND_MAX;
      //a[i * p + j] = (double) (i+j);
      //printf("%.2f  ", a[i * p + j]);
    }
    //printf("\n");
  }

  //printf("b\n");
  for(i=0;i < p;i++){
    for(j=0;j < m;j++) {
      b[i * m + j] = (double) rand() / (double) RAND_MAX;
      //b[i * m + j] = (double) (i+j);
      //printf("%.2f  ", b[i * m + j]);
    }
    //printf("\n");
  }

  // ------------- COMPUTATION DONE ON GPU ----------------------------

  cudaMalloc((void**)&dev_a, size_a); // allocate memory on device
  cudaMalloc((void**)&dev_b, size_b);
  cudaMalloc((void**)&dev_c, size_c);

  cudaMemcpy(dev_a, a , size_a ,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b , size_b ,cudaMemcpyHostToDevice);

  cudaEventCreate(&start); // instrument code to measure start time
  cudaEventCreate(&stop);
  
  cudaEventRecord(start, 0);
  // cudaEventSynchronize(start); // not needed

  gpu_matrixmult<<<Grid,Block>>>(dev_a,dev_b,dev_c,n,m,p);

  cudaEventRecord(stop, 0); // instrument code to measure end time
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time_ms, start, stop );

  cudaMemcpy(c,dev_c, size_c ,cudaMemcpyDeviceToHost);

  printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time
  /*
  printf("c\n");
  for(i=0;i < n;i++){
    for(j=0;j < m;j++) {
      printf("%.2f  ", c[i * m + j]);
    }
    printf("\n");
  }*/
  
  // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
  // DEBUGGING USE ONLY (AND FOR LIMITED NUMBERS OF TIMING RUNS)
/*
  cudaEventRecord(start, 0); // use same timing
  // cudaEventSynchronize(start); // not needed

  cpu_matrixmult(a,b,c, n, m, p); // do calculation on host (NOTE: This computes the diff with GPU result.)

  cudaEventRecord(stop, 0); // instrument code to measue end time
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time_ms, start, stop );

  printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms); // exec. time

// ------------------- check device creates correct results -----------------

  double error, suma, sumb, sumc, ai, bi, ci;
  suma = 0.; sumb = 0; sumc = 0;
  for(i=0;i < n*n;i++) {
    ai = (double) a[i];
    bi = (double) b[i];
    ci = (double) c[i];
    suma += ai*ai;
    sumb += bi*bi;
    sumc += ci*ci;
  }
  suma = sqrt(suma);
  sumb = sqrt(sumb);
  sumc = sqrt(sumc);
  error =  sumc/(n*suma*sumb);
  printf("Scaled error between GPU and CPU: %e\n", error);
  
*/
// -------------- clean up ---------------------------------------

  free(a);
  free(b);
  free(c);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  return 0;
}
