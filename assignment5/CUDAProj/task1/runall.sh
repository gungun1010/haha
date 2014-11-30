#!/bin/bash

./Task1CPUsp 1024 1024 1024 > outputs/Task1CPUsp_1024_1024_1024.out
printf "1 "
./Task1CPUsp 2048 2048 2048 > outputs/Task1CPUsp_2048_2048_2048.out
printf "2 "
./Task1CPUsp 4096 4096 4096 > outputs/Task1CPUsp_4096_4096_4096.out
printf "3 "
./Task1CPUsp 8192 8192 8192 > outputs/Task1CPUsp_8192_8192_8192.out
printf "4 "
./Task1CPUsp 1024 1024 8192 > outputs/Task1CPUsp_1024_1024_8192.out
printf "5 "
./Task1CPUsp 8192 8192 1024 > outputs/Task1CPUsp_8192_8192_1024.out
printf "6 "
./Task1CPUsp 8192 1024 8192 > outputs/Task1CPUsp_8192_1024_8192.out
printf "7 "
./Task1CPUdp 8192 8192 8192 > outputs/Task1CPUdp_1024_1024_1024.out
