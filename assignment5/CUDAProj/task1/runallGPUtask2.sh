#!/bin/bash
./Task2GPUsp 1024 1024 1024 32 32 32 32 32 > outputs/Task2GPUsp_1024_1024_1024.out
printf "1 "
./Task2GPUsp 2048 2048 2048 32 32 64 64 32 > outputs/Task2GPUsp_2048_2048_2048.out
printf "2 "
./Task2GPUsp 4096 4096 4096 32 32 128 128 32 > outputs/Task2GPUsp_4096_4096_4096.out
printf "3 "
./Task2GPUsp 8192 8192 8192 32 32 256 256 32 > outputs/Task2GPUsp_8192_8192_8192.out
printf "4 "
./Task2GPUsp 1024 1024 8192 32 32 32 32 32 > outputs/Task2GPUsp_1024_1024_8192.out
printf "5 "
./Task2GPUsp 8192 8192 1024 32 32 256 256 32 > outputs/Task2GPUsp_8192_8192_1024.out
printf "6 "
./Task2GPUsp 8192 1024 8192 32 32 32 256 32 > outputs/Task2GPUsp_8192_1024_8192.out
printf "7 "
./Task3GPUsp 8192 8192 8192 32 32 256 256 32 > outputs/Task3GPUsp_8192_8192_8192.out
