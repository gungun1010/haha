#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>
#include <math.h>

__global__ void matmul_tile(float *a, float *b, float *c, int n, int m, int p, int TW, int NTB) {
  extern __shared__ float bigarray[]; 

  float *aTile=&bigarray[0], *bTile=&bigarray[TW*TW];
  int tx = threadIdx.x; 
  int ty = threadIdx.y; 
  float *cvalue;//scope: thread 
  int col = tx + blockDim.x * blockIdx.x;
  int row = ty + blockDim.y * blockIdx.y;
  int tileNum, aIdx, bIdx, tileIdx_m, tileCol;

  tileNum = p/TW + (p % TW != 0);

  cvalue = (float *)malloc(NTB*sizeof(float));

  //init c tiles
  for (tileIdx_m=0; tileIdx_m<NTB; tileIdx_m++) cvalue[tileIdx_m] = 0.;

  for (int tileIdx_p=0; tileIdx_p<tileNum; tileIdx_p++) {
    
    //load aTile
    aIdx = tileIdx_p*TW + tx;
    if(aIdx >= p || row >= n){
        aTile[ty*TW+tx] = 0.;
    }else{
        aTile[ty*TW+tx] = a[row*p + aIdx]; //Copy to shared memory 
    }
    
    for(tileIdx_m=0; tileIdx_m<NTB; tileIdx_m++){ 
        //load btile[ty][tx] with element [ty][tx] in tileIdx_m-th tile of b
        if((blockIdx.x % NTB) == 0){
            bIdx = tileIdx_p*TW +ty;
            tileCol = tx + blockDim.x*(blockIdx.x+tileIdx_m);
            if(bIdx >= p || tileCol >= m){
                bTile[ty*TW+tx] = 0.;
            }else{
                bTile[ty*TW+tx] = b[bIdx*m + tileCol]; //Copy to shared memory 
            }

            __syncthreads();
            for (int k=0; k<TW; k++){
                 cvalue[tileIdx_m] += aTile[ty*TW+k] * bTile[k*TW+tx];
                 //printf("bx = %d, by = %d, (tx = %d, ty = %d) @ tileIdx_m = %d : a=%.2f b=%.2f \n",blockIdx.x, blockIdx.y, tx, ty, tileIdx_m, aTile[ty*TW+k],bTile[k*TW+tx]);
            }
            //printf("bx = %d, by = %d, (tx = %d, ty = %d) @ tileIdx_m = %d: c= %.2f\n",blockIdx.x, blockIdx.y, tx, ty, tileIdx_m, cvalue[tileIdx_m]);
            __syncthreads();
            c[row*m + tileCol] = cvalue[tileIdx_m];
        }

    }
  }

  if(row < n && col < m){
    for(tileIdx_m=0; tileIdx_m<NTB; tileIdx_m++){ 
        //load to C
        if((blockIdx.x % NTB) == 0){
            tileCol = tx + blockDim.x*(blockIdx.x+tileIdx_m);
            c[row*m + tileCol] = cvalue[tileIdx_m];
        }
    }
  }
    
  free(cvalue);
}


void cpu_matrixmult(float *a,float *b, float *c, int n, int m, int p) {

  int index, indexa, indexb;
  float cvalue;
  for(int col=0;col < m; col++){
    for(int row=0;row < n; row++) {
      indexb = col;
      index = row * m + col;
      cvalue = 0.;
      for (indexa = row*p; indexa < (row*p + p); indexa++, indexb+=m){ 
        cvalue += a[indexa]*b[indexb];
      }
      c[index] -= cvalue; //NOTE: This calculates the diff between CPU and GPU computations.
    }
  }
}


int main(int argc, char *argv[]) {

  int i, j; // loop counters

  int gpucount = 0; // Count of available GPUs
  int Grid_Dim_x = 1; //Grid dimension, x
  int Grid_Dim_y = 1; //Grid dimension, y
  int Block_Dim_x = 1; //Block dimension, x
  int Block_Dim_y = 1; //Block dimension, y
  int TW = 1;
  int NTB = 1;

  int n,m,p; // matrix dimension
  float *a,*b,*c;
  float *dev_a, *dev_b, *dev_c;
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
  if (argc<10) {
    printf("# of inputs: %d\n", argc);
    printf("Usage: Task1GPUsp <n> <m> <p> <block dim x> <block dim y> <grid dim x> <grid dim y> <tile width> <Number of tiles>\n");
    exit (-1);
  }

  n = atoi(argv[1]);
  m = atoi(argv[2]);
  p = atoi(argv[3]);
  

  Block_Dim_x = atoi(argv[4]); // non-Square block, # of rows
  Block_Dim_y = atoi(argv[5]); // non-Square block, # of cols
  if (Block_Dim_x * Block_Dim_y > 1024) {
    printf("Error, too many threads in block\n");
    exit (-1);
  }

  //not really used in Task2 
  Grid_Dim_x = atoi(argv[6]); // non-Square grid, # of rows
  Grid_Dim_y = atoi(argv[7]); // non-Square grid, # of cols
  
  TW = atoi(argv[8]);
   
  if(Block_Dim_x != Block_Dim_y || Block_Dim_x != TW || Block_Dim_y != TW){
      printf("Error, bx, by, tw must be equal\n");
      exit(-1);
  }

  //printf("A Matrix Dimension = %dx%d\n",n,p);
  //printf("B Matrix Dimension = %dx%d\n",p,m);
  //printf("C Matrix Dimension = %dx%d\n",n,m);
  Grid_Dim_x = m/Block_Dim_x + (m % Block_Dim_x != 0);
  Grid_Dim_y = n/Block_Dim_y + (n % Block_Dim_y != 0);

  NTB = atoi(argv[9]);

  //printf("Grid_x = %d Grid_y = %d NTB = %d\n", Grid_Dim_x,Grid_Dim_y,NTB);

  dim3 Grid(Grid_Dim_x, Grid_Dim_y); //Grid structure
  dim3 Block(Block_Dim_x, Block_Dim_y); //Block structure

  size_a = n * p * sizeof(float); // number of bytes in total in arrays
  size_b = p * m * sizeof(float); // number of bytes in total in arrays
  size_c = n * m * sizeof(float); // number of bytes in total in arrays

  a = (float*) malloc(size_a); // dynamically allocated memory for arrays on host
  b = (float*) malloc(size_b);
  c = (float*) malloc(size_c); // results from GPU

  srand(12345);
  //int p = n; //Used here only to illustrate proper initialization for non-square case
   
  //printf ("a\n");
  for(i=0;i < n;i++){
    for(j=0;j < p;j++) {
      a[i * p + j] = (float) rand() / (float) RAND_MAX;
      //a[i * p + j] = (float) (i+j);
      //printf("%.2f  ", a[i * p + j]);
    }
    //printf("\n");
  }

  //printf("b\n");
  for(i=0;i < p;i++){
    for(j=0;j < m;j++) {
      b[i * m + j] = (float) rand() / (float) RAND_MAX;
      //b[i * m + j] = (float) (i+j);
      //printf("%.2f  ", b[i * m + j]);
    }
    //printf("\n");
  }

  // ------------- COMPUTATION DONE ON GPU ----------------------------

  errorcode = cudaMalloc((void**)&dev_a, size_a); // allocate memory on device
    if(errorcode != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("cudaMalloc error: %s\n", cudaGetErrorString(errorcode));
        exit(-1);
    }
  errorcode = cudaMalloc((void**)&dev_b, size_b);
    if(errorcode != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("cudaMalloc error: %s\n", cudaGetErrorString(errorcode));
        exit(-1);
    }
  errorcode = cudaMalloc((void**)&dev_c, size_c);
    if(errorcode != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("cudaMalloc error: %s\n", cudaGetErrorString(errorcode));
        exit(-1);
    }

  cudaMemcpy(dev_a, a , size_a ,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b , size_b ,cudaMemcpyHostToDevice);

  cudaEventCreate(&start); // instrument code to measure start time
  cudaEventCreate(&stop);
  
  cudaEventRecord(start, 0);
  // cudaEventSynchronize(start); // not needed
  size_t Ns = 2 * TW*TW * sizeof(float);
  size_t heapSize = Grid_Dim_x * Grid_Dim_y * Block_Dim_x* Block_Dim_y * NTB * sizeof(float)/4; 
  errorcode = cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapSize);
    if(errorcode != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("cuda device heap error: %s\n", cudaGetErrorString(errorcode));
        exit(-1);
    }

  matmul_tile<<<Grid,Block, Ns>>>(dev_a, dev_b, dev_c, n, m, p, TW, NTB);

    // make the host block until the device is finished with foo
    cudaThreadSynchronize();

    // check for error
    errorcode = cudaGetLastError();
    if(errorcode != cudaSuccess)
    {
        // print the CUDA error message and exit
        printf("CUDA error: %s\n", cudaGetErrorString(errorcode));
        exit(-1);
    }

  cudaEventRecord(stop, 0); // instrument code to measure end time
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time_ms, start, stop );

  cudaMemcpy(c,dev_c, size_c ,cudaMemcpyDeviceToHost);

  //printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms); // exec. time
   
  //printf("c\n");
  for(i=0;i < n;i++){
    for(j=0;j < m;j++) {
      printf("%.2f  ", c[i * m + j]);
    }
    printf("\n");
  }
  
  // ------------- COMPUTATION DONE ON HOST CPU ----------------------------
  // DEBUGGING USE ONLY (AND FOR LIMITED NUMBERS OF TIMING RUNS)
/*
  cudaEventRecord(start, 0); // use same timing

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
