#!/bin/bash
icc -o pi_o0 -g -O0 -fno-alias -std=c99 pi.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o pi_o1 -g -O1 -fno-alias -std=c99 pi.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o pi_o3_raw -g -O3 -no-vec -no-simd -fno-alias -std=c99 pi.c ~ahs3/cpsc424/utils/timing/timing.o
icc -o pi_o3 -g -O3 -fno-alias -std=c99 pi.c ~ahs3/cpsc424/utils/timing/timing.o
