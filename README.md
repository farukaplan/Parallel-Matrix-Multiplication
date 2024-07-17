# Parallel Matrix Multiplication with MPI

## Build and Run
```sh	
# Sequential program
mpicc -o sequential sequential.c
mpiexec -n 1 ./sequential data/BigA.txt data/BigX.txt
```

```sh
# Parallel p2p program
mpicc -o p2p point-to-point.c
mpiexec -n 4 ./p2p data/BigA.txt data/BigX.txt
```

```sh
# Parallel collvetive program
mpicc -o collective collective.c
mpiexec -n 4 ./collective data/BigA.txt data/BigX.txt
```

You can apply this to desired matrices, by simply changing the input txt names
- Don't forget, the matrix dimensions should be match, otherwise it does not work

Also you can adjust the core number for parallel program as you wish
- The core number can be 1, 2, 4, 8, 16 and goes on
- Don't exceed your computer's limitations

## Discussion
You can find an explanation of codes and discussion about how execution time changes when we change core number, and other stuff in Report.pdf

## To Do
- [ ] Write makefile
