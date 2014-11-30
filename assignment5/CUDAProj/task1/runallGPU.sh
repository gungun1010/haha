#!/bin/bash
./Task1GPUsp 1024 1024 1024 32 32 32 32 > outputs/Task1GPUsp_1024_1024_1024.out
printf "1 "
./Task1GPUsp 2048 2048 2048 32 32 64 64 > outputs/Task1GPUsp_2048_2048_2048.out
printf "2 "
./Task1GPUsp 4096 4096 4096 32 32 128 128 > outputs/Task1GPUsp_4096_4096_4096.out
printf "3 "
./Task1GPUsp 8192 8192 8192 32 32 256 256 > outputs/Task1GPUsp_8192_8192_8192.out
printf "4 "
./Task1GPUsp 1024 1024 8192 32 32 32 32 > outputs/Task1GPUsp_1024_1024_8192.out
printf "5 "
./Task1GPUsp 8192 8192 1024 32 32 256 256 > outputs/Task1GPUsp_8192_8192_1024.out
printf "6 "
./Task1GPUsp 8192 1024 8192 32 32 32 256  > outputs/Task1GPUsp_8192_1024_8192.out
printf "7 "
./Task1GPUdp 8192 8192 8192 32 32 256 256 > outputs/Task1GPUdp_1024_1024_1024.out
