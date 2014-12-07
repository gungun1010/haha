/*
 *  Compilation: gcc -Wall ex1.c -o ex1 -L/home/leon/CUDAProj/finalProj/clamav/lib -lclamav
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include </home/leon/clamav/include/clamav.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define TAG_BITS 10

const char *BYTE_CODE = "bytecode.cvd";
const char *DAILY = "/dailyPack/dailyGPUsig.bin";
const char *MAIN = "/mainPack/mainGPUsig.bin";

__global__ void patternMatching(char *set1, char *set2, char *set3){
}

//function to construct the signature database for GPU
void loadSig (const char *fileName, char **buffer, size_t *size){
    long lSize;
    FILE *fp;
    char *fullPath; 
    
    fullPath = (char *)malloc(100*sizeof(char));

    //find file
    strcpy(fullPath, (char *) cl_retdbdir());//this function returns where the signature file is
    strcat(fullPath, "/");
    strcat(fullPath, fileName);
    fp = fopen (fullPath , "rb" );
    if( !fp ) 
        perror(fileName),exit(1);

    //seek the end of file
    fseek( fp , 0L , SEEK_END);
    lSize = ftell( fp );
    rewind( fp );
    //printf("%ld\n",lSize);
    (*size) = lSize;

    /* allocate memory for entire content */
    (*buffer) = (char *) calloc( 1, lSize+1 );
    if( !(*buffer) ) 
        fclose(fp),fputs("memory alloc fails",stderr),exit(1);

    /* copy the file into the buffer */
    if( 1!=fread( (*buffer) , lSize, 1 , fp) )
        fclose(fp),free((*buffer)),fputs("entire read fails",stderr),exit(1);

    fclose(fp);
    //free(buffer);
}

//function to load input file to scan
void loadFile (const char *fileName, char **buffer, size_t *size){
    long lSize;
    FILE *fp;
    char *fullPath; 
    
    fullPath = (char *)malloc(100*sizeof(char));
    //find file
    strcpy(fullPath, (char *) cl_retdbdir());//this function returns where the signature file is
    //printf("%s\n",fullPath);
    strcat(fullPath, "/");
    //printf("%s\n",fullPath);
    strcat(fullPath, fileName);
    //printf("find signature file in %s\n",fullPath);
    fp = fopen (fullPath , "rb" );
    if( !fp ) perror(fileName),exit(1);
    
    //seek the beginning of file
    //fseek(fp, SEEK_SET, 0);
    fseek( fp , 0L , SEEK_END);
    lSize = ftell( fp );
    rewind( fp );
    //printf("%ld\n",lSize);
    (*size) = lSize;

    /* allocate memory for entire content */
    (*buffer) = (char *) calloc( 1, lSize+1 );
    if( !(*buffer) ) 
        fclose(fp),fputs("memory alloc fails",stderr),exit(1);

    /* copy the file into the buffer */
    if( 1!=fread( (*buffer) , lSize, 1 , fp) )
          fclose(fp),free((*buffer)),fputs("entire read fails",stderr),exit(1);

    fclose(fp);
    //free(buffer);
}
/*
 * Exit codes:
 *  0: clean
 *  1: infected
 *  2: error
 */
//const char *DBDIR = "/home/leon/clamav/share/clamav";

int main(int argc, char **argv)
{
	int fd, ret;
	unsigned long int size = 0;
	unsigned int sigs = 0;
	long double mb;
	const char *virname;
	struct cl_engine *engine;

    int gpucount = 0; // Count of available GPUs
    //We only have 3701312 signatures
    //each thread get 1 signature, we need no more than 1024*1024 threads
    //grid size is then fixed to (32,32,1), and block size is (32,32,1)

    int Grid_Dim_x = 256; //Grid dimension, x
    int Grid_Dim_y = 256; //Grid dimension, y
    int Block_Dim_x = 32; //Block dimension, x
    int Block_Dim_y = 32; //Block dimension, y
    cudaError_t errorcode;

    //host buffer to store each signature dataset
    char *byteCodeBuf; 
    char *dailyBuf;
    char *mainBuf;
    char *devBcb, *devDb, *devMb;//device buffer correspoding to the host buffer
    size_t sizeBcb, sizeDb, sizeMb;
   
    // --------------------SET PARAMETERS AND DATA -----------------------
    //load signatures into host buffer
    //loadFile(BYTE_CODE, &byteCodeBuf, &sizeBcb);
    loadSig(DAILY, &dailyBuf, &sizeDb);
    loadSig(MAIN, &mainBuf, &sizeMb);
    
    for(int i=0; i<11; i++){
        printf("%x ", (unsigned char) dailyBuf[i]);
    }

    exit(1);
    //loadFile(MAIN, &mainBuf, &sizeMb);

    errorcode = cudaGetDeviceCount(&gpucount);
    if (errorcode == cudaErrorNoDevice) {
        printf("No GPUs are visible\n");
        exit(-1);
    }

    if (Block_Dim_x * Block_Dim_y > 1024) {
        printf("Error, too many threads in block\n");
        exit (-1);
    }

    dim3 Grid(Grid_Dim_x, Grid_Dim_y); //Grid structure
    dim3 Block(Block_Dim_x, Block_Dim_y); //Block structure
    
    cudaMalloc((void**)&devBcb, sizeBcb*sizeof(char));
    cudaMalloc((void**)&devDb, sizeDb*sizeof(char));
    cudaMalloc((void**)&devMb, sizeMb*sizeof(char)); 

    cudaMemcpy(devBcb, byteCodeBuf , sizeBcb ,cudaMemcpyHostToDevice);
    cudaMemcpy(devDb, dailyBuf , sizeDb ,cudaMemcpyHostToDevice);
    cudaMemcpy(devMb, mainBuf , sizeMb ,cudaMemcpyHostToDevice);

    if(argc != 2) {
        printf("Usage: %s file\n", argv[0]);
        return 2;
    }

    if((fd = open(argv[1], O_RDONLY)) == -1) {
        printf("Can't open file %s\n", argv[1]);
        return 2;
    }


    if((ret = cl_init(CL_INIT_DEFAULT)) != CL_SUCCESS) {
        printf("Can't initialize libclamav: %s\n", cl_strerror(ret));
        return 2;
    }

    if(!(engine = cl_engine_new())) {
        printf("Can't create new engine\n");
        return 2;
    }
    /* load all available databases from default directory */
    printf("loading signatures in %s\n",cl_retdbdir());
    if((ret = cl_load(cl_retdbdir(), engine, &sigs, CL_DB_STDOPT)) != CL_SUCCESS) {
        printf("cl_load: %s\n", cl_strerror(ret));
        close(fd);
            cl_engine_free(engine);
        return 2;
    }

    printf("Loaded %u signatures.\n", sigs);

    /* build engine */
    if((ret = cl_engine_compile(engine)) != CL_SUCCESS) {
        printf("Database initialization error: %s\n", cl_strerror(ret));;
            cl_engine_free(engine);
        close(fd);
        return 2;
    }

    /* scan file descriptor */
    if((ret =cl_scandesc(fd, &virname, &size, engine, CL_SCAN_STDOPT)) == CL_VIRUS) {
        printf("Virus detected: %s\n", virname);
    } else {
        if(ret == CL_CLEAN) {
            printf("No virus detected.\n");
        } else {
            printf("Error: %s\n", cl_strerror(ret));
            cl_engine_free(engine);
            close(fd);
            return 2;
        }
    }
    close(fd);

    /* free memory */
    cl_engine_free(engine);

    /* calculate size of scanned data */
    mb = size * (CL_COUNT_PRECISION / 1024) / 1024.0;
    printf("Data scanned:%ld  %2.2Lf MB\n", size, mb);

    return ret == CL_VIRUS ? 1 : 0;
}
