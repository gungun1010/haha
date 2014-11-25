CUDAPATH = /usr/local/cuda-6.5

NVCC = $(CUDAPATH)/bin/nvcc

NVCCFLAGS = -I$(CUDAPATH)/include -O3

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

# Compiler-specific flags (by default, we always use sm_20)
GENCODE_SM20 = -gencode=arch=compute_50,code=\"sm_50,compute_50\"
GENCODE = $(GENCODE_SM20)

.SUFFIXES : .cu .ptx .c

BINARIES = Task1GPUsp Task1CPUsp

all: Task1GPUsp Task1CPUsp

Task1GPUsp: Task1GPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task1CPUsp: Task1CPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

.cu.o:
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

clean:	
	rm -f *.o $(BINARIES)
