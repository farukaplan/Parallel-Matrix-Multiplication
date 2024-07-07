// Faruk KAPLAN
// 21050111026

// For compile: mpicc -o p2p point-to-point.c
// For run: mpiexec -n 4 ./p2p BigA.txt BigX.txt (you can change core number or file names)

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

double* multiply(double *portion_A, double **x, int size, int col);
double* flatten(double** matrix, int row, int col);
void readMatrix(const char *filename, double ***matrix, int *row, int *col);
void freeMatrix(double **matrix, int rows);

int main(int argc, char** argv) {
    if(argc != 3) {
        printf("A few or a lot argument !");
        return 1;
    }

    // Variables for measure time
    clock_t local_start, local_end; 
    double local_total, total_time;

    // Take file names as command-line argument
    const char *filename_A = argv[1];
    const char *filename_x = argv[2];

    // Our matrices variables
    double *result_portion, *result_portion_A, *result_portion_process;
    int row_A, col_A, row_x, col_x;
    double **A, **x;

    // Read the content of the files into 2d arrays
    readMatrix(filename_A, &A, &row_A, &col_A);
    readMatrix(filename_x, &x, &row_x, &col_x);

    // Flatten A matrix to send and receive
    double *flattened_A = flatten(A, row_A, col_A);

    int proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // We want to split the matrix equally 
    if (row_A % proc != 0) {
        if (rank == 0) {
            printf("Matrix size must be divisible by the number of processes.\n");
        }

        MPI_Finalize();
        return 1;
    }

    // Master core
    if(rank == 0) {
        // Size of row for each process
        int row_process = row_A / proc;

        // Allocate space for temporary matrix
        double *tempMatrix = (double *)malloc(row_process * col_A * sizeof(double));

        // Send each process to its portion
        for(int i = 1; i < proc; i++) {
            for(int j = 0; j < (row_process * col_A); j++) {
                tempMatrix[j] = flattened_A[(i * row_process * col_A) + j];
            }

            MPI_Send(tempMatrix, row_process * col_A, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        // This process will also do computation with its own portion
        for(int i = 0; i < row_process * col_A; i++) {
            tempMatrix[i] = flattened_A[i];
        }

        result_portion_A = (double *)malloc(row_process * sizeof(double));

        // This process's result, also we'll measure the time
        local_start = clock();
        result_portion_A = multiply(tempMatrix, x, row_process, col_A);
        local_end = clock();

        // Only master processes time
        total_time = ((double) (local_end - local_start)) / CLOCKS_PER_SEC;

        // Printing the result onto file
        FILE *output = fopen("y.txt", "w");

        for(int i = 0; i < row_process; i++) {
            fprintf(output, "%.2lf\n", result_portion_A[i]);
        }

        // An array for collectiong results of other processes
        result_portion_process = (double *)malloc(row_process * sizeof(double));

        // Collect other processes results
        for(int i = 1; i < proc; i++) {
            MPI_Recv(result_portion_process, row_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for(int j = 0; j < row_process; j++) {
                fprintf(output, "%.2lf\n", result_portion_process[j]);
            }

            // Receive the time
            double time = 0.0;
            MPI_Recv(&time, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_time += time;
        }

        printf("Time taken for multiplication in p2p mode: %lf seconds\n", total_time);
        
    // Other processes
    } else {
        int row_process = row_A / proc;

        // Array for receiving the portion
        double *portion_A = (double *)malloc(row_process * col_A * sizeof(double));

        // Each process receives its own portion
        MPI_Recv(portion_A, row_process * col_A, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        result_portion = (double *)malloc(row_process * sizeof(double));

        // Each process calculates its own results, also we'll measure the local time
        local_start = clock();
        result_portion = multiply(portion_A, x, row_process, col_A);
        local_end = clock();

        // Calculate total local time, sum them and send it to master process
        local_total = ((double) (local_end - local_start)) / CLOCKS_PER_SEC;

        // Each process sends their results to master process
        MPI_Send(result_portion, row_process, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        // Send the time value
        MPI_Send(&local_total, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    // Freeing up A and x matrices
    freeMatrix(A, row_A);
    freeMatrix(x, row_x);

    return 0;
}

// Multiplication
double* multiply(double *portion_A, double **x, int row, int col) {
    double* result = (double *)malloc(row * sizeof(double));

    for(int i = 0; i < row; i++) {
        double tempResult = 0.0;

        for(int j = 0; j < col; j++) {
            tempResult += portion_A[(i * col) + j] * x[j][0];
        }

        result[i] = tempResult;
    }

    return result;
}

// Flatten the matrix
double* flatten(double** matrix, int row, int col) {
    double* flattened = (double*)malloc(row * col * sizeof(double));
    int k = 0;

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            flattened[k++] = matrix[i][j];
        }
    }

    return flattened;
}

// Reading the matrix from file
void readMatrix(const char *filename, double ***matrix, int *row, int *col) {
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    // Read the dimensions
    fscanf(file, "%d %d", row, col);

    // Allocate memory for matrix
    *matrix = (double **)malloc((*row) * sizeof(double *));

    if (*matrix == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < *row; i++) {
        (*matrix)[i] = (double *)malloc((*col) * sizeof(double));
        if ((*matrix)[i] == NULL) {
            printf("Memory allocation failed\n");
            exit(1);
        }
    }

    for (int i = 0; i < *row; i++) {
        for (int j = 0; j < *col; j++) {
            fscanf(file, "%lf", &(*matrix)[i][j]);
        }
    }

    fclose(file);
}

// Free allocated memory
void freeMatrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}