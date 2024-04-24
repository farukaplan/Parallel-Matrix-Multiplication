# Parallel-Matrix-Multiplication

This repo contains 3 main files: 
1-) Sequential.c: This c code multiplies matrices in traditional, sequential way
2-) Point-to-point.c: This c code is parallelized version of first file, but only using point to point communication protocols (MPI_Send and MPI_Recv)
3-) Collective.c: This c code is also parallelized version of first file, but also using collective MPI functions

There is an information at the beginning of the each file that how you can compile and run files

When you run the code, you can change input file names and core number, input file format should be same as my input files, also give proper number for cores
