// Faruk KAPLAN
// 21050111026

// For compile: mpicc -o sequential sequential.c
// For run: mpiexec -n 1 ./sequential BigA.txt BigX.txt (you can change the file names)

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

void multiply(double **A, double **x, double **result, int row_A, int col_A, int row_x, int col_x);
void readMatrix(const char *filename, double ***matrix, int *row, int *col);
void freeMatrix(double **matrix, int rows);

int main(int argc, char** argv) {
    if(argc != 3) {
        printf("A few or a lot argument !");
        return 1;
    }
    
    // We want to compile with mpicc, not gcc (for correct comparison)
    int proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Take file names as command-line argument
    const char *filename_A = argv[1];
    const char *filename_x = argv[2];

    // Our matrices variables
    int row_A, col_A, row_x, col_x;
    double **A, **x;

    // Read the content of the files into 2d arrays
    readMatrix(filename_A, &A, &row_A, &col_A);
    readMatrix(filename_x, &x, &row_x, &col_x);

    // Allocate space for result matrix
    double **result = (double **)malloc(row_A * sizeof(double *)); 
    for (int i = 0; i < row_A; i++) {
        (result)[i] = (double *)malloc(col_x * sizeof(double));
    }
    
    // Start to measure time
    clock_t start = clock();
    multiply(A, x, result, row_A, col_A, row_x, col_x);
    clock_t end = clock();

    double time_spent = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for multiplication in sequential mode: %lf seconds\n", time_spent);

    // Writing output matrix into output file
    FILE *output = fopen("y.txt", "w");

    for (int i = 0; i < row_A; i++) {
        for (int j = 0; j < col_x; j++) {
            fprintf(output, "%.2lf\t", result[i][j]);
        }
        fprintf(output, "\n");
    }

    // Freeing up allocated matrices
    freeMatrix(A, row_A);
    freeMatrix(x, row_x);

    MPI_Finalize();
    
    return 0;
}

// Perform multiplication
void multiply(double **A, double **x, double **result, int row_A, int col_A, int row_x, int col_x) {
    int i, j, k;

    for (i = 0; i < row_A; i++) {
        for (j = 0; j < col_x; j++) {
            result[i][j] = 0.0;

            for (k = 0; k < col_A; k++) {
                result[i][j] += A[i][k] * x[k][j];
            }
        }
    }
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