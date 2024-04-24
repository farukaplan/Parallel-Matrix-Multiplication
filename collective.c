// Faruk KAPLAN
// 21050111026

// For compile: mpicc -o collective collective.c
// For run: mpiexec -n 4 ./collective BigA.txt BigX.txt (you can change core number or file names)

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

    // Take file names as command-line argument
    const char *filename_A = argv[1];
    const char *filename_x = argv[2];

    // Variables for measure time
    clock_t local_start, local_end;

    // Our matrices variables
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

    // Calculate number of rows per process
    int row_process = row_A / proc;

    // Allocate memory for the portion of A
    double *portion_A = (double*)malloc(row_process * col_A * sizeof(double));

    // Scatter flattened A to all processes
    MPI_Scatter(flattened_A, row_process * col_A, MPI_DOUBLE, portion_A, row_process * col_A, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Perform multiplication and measure time
    local_start = clock();
    double *result = multiply(portion_A, x, row_process, col_A);
    local_end = clock();

    double total_time = 0.0;
    total_time += ((double)(local_end - local_start)) / CLOCKS_PER_SEC;

    // Gather results back to the master process
    double *gathered_results = NULL;

    // Only allocate in master process
    if (rank == 0) {
        gathered_results = (double*)malloc(row_A * sizeof(double));
    }

    // Gather the results
    MPI_Gather(result, row_process, MPI_DOUBLE, gathered_results, row_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Output or process gathered results as needed
    if (rank == 0) {
        FILE *output = fopen("y.txt", "w");

        for(int i = 0; i < row_A; i++) {
            fprintf(output, "%.2lf\n", gathered_results[i]);
        }

        printf("Time taken for multiplication in collective mode: %lf seconds\n", total_time);
    }

    // Free allocated memory
    free(portion_A);
    free(result);
    
    if (rank == 0) {
        free(gathered_results);
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