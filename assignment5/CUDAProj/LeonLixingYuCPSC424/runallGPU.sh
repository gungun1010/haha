#!/bin/bash
make gpu #Make sure your binaries have the same name as shown here
./Task1GPUsp 1024 1024 1024 32 32 32 32 > Task1GPUsp.out
./Task1GPUdp 1024 1024 1024 32 32 32 32 > Task1GPUdp.out
./Task2GPUsp 1024 1024 1024 32 32 32 32 32 > Task2GPUsp.out
./Task2GPUdp 1024 1024 1024 32 32 32 32 32 > Task2GPUdp.out
./Task3GPUsp 1024 1024 1024 32 32 32 32 32 32 > Task3GPUsp.out
