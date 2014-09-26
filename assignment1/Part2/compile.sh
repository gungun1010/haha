#!/bin/bash
icc -o vectorTriad_o0 -g -O0 -fno-alias -std=c99 vectorTriad.c dummy.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o vectorTriad_o1 -g -O1 -fno-alias -std=c99 vectorTriad.c dummy.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o vectorTriad_o3_raw -g -O3 -no-vec -no-simd -fno-alias -std=c99 vectorTriad.c dummy.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o vectorTriad_o3 -g -O3 -xHost -fno-alias -std=c99 vectorTriad.c dummy.c ~ahs3/cpsc424/utils/timing/timing.o
