CUDAPATH = /home/apps/fas/GPU/cuda_6.0.37

NVCC = $(CUDAPATH)/bin/nvcc

NVCCFLAGS = -I$(CUDAPATH)/include -O3

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

# Compiler-specific flags (by default, we always use sm_20)
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

.SUFFIXES : .cu .ptx .c

BINARIES = Task1GPUsp Task1CPUsp Task1GPUdp Task1CPUdp Task2GPUsp Task2GPUdp Task3GPUsp

all: $(BINARIES)

cpu: Task1CPUsp Task1CPUdp
gpu: Task1GPUsp Task1GPUdp Task2GPUsp Task2GPUdp Task3GPUsp


Task1GPUsp: Task1GPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task1CPUsp: Task1CPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task1GPUdp: Task1GPUdp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task1CPUdp: Task1CPUdp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task2GPUsp: Task2GPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task2GPUdp: Task2GPUdp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

Task3GPUsp: Task3GPUsp.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<

.cu.o:
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

clean:	
	rm -f *.o *.out $(BINARIES)
