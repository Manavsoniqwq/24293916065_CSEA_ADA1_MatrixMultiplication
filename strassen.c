#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// --- Matrix Memory and Utility ---

int **create_square_matrix(int dimension_n) {
    int **matrix_handle = (int **)malloc(dimension_n * sizeof(int *));
    for (int i = 0; i < dimension_n; i++) {
        matrix_handle[i] = (int *)malloc(dimension_n * sizeof(int));
    }
    return matrix_handle;
}

void destroy_matrix(int dimension_n, int **matrix_handle) {
    for (int i = 0; i < dimension_n; i++) {
        free(matrix_handle[i]);
    }
    free(matrix_handle);
}

void display_matrix(int dimension_n, int **matrix) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void add_matrices(int dimension_n, int **MatrixA, int **MatrixB, int **MatrixResult) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            MatrixResult[i][j] = MatrixA[i][j] + MatrixB[i][j];
        }
    }
}

void subtract_matrices(int dimension_n, int **MatrixA, int **MatrixB, int **MatrixResult) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            MatrixResult[i][j] = MatrixA[i][j] - MatrixB[i][j];
        }
    }
}

// --- Strassen's O(n^2.807) Algorithm ---

void multiply_strassen_recursive(int n, int **A, int **B, int **C) {
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    
    int sub_n = n / 2;
    
    // Allocate 19 temporary matrices
    int **A11 = create_square_matrix(sub_n); int **A12 = create_square_matrix(sub_n);
    int **A21 = create_square_matrix(sub_n); int **A22 = create_square_matrix(sub_n);
    int **B11 = create_square_matrix(sub_n); int **B12 = create_square_matrix(sub_n);
    int **B21 = create_square_matrix(sub_n); int **B22 = create_square_matrix(sub_n);
    int **C11_res = create_square_matrix(sub_n); int **C12_res = create_square_matrix(sub_n);
    int **C21_res = create_square_matrix(sub_n); int **C22_res = create_square_matrix(sub_n);
    int **P1 = create_square_matrix(sub_n); int **P2 = create_square_matrix(sub_n);
    int **P3 = create_square_matrix(sub_n); int **P4 = create_square_matrix(sub_n);
    int **P5 = create_square_matrix(sub_n); int **P6 = create_square_matrix(sub_n);
    int **P7 = create_square_matrix(sub_n);
    int **TempMatrix1 = create_square_matrix(sub_n); 
    int **TempMatrix2 = create_square_matrix(sub_n);
    
    // Partition A and B
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            A11[i][j] = A[i][j]; A12[i][j] = A[i][j + sub_n];
            A21[i][j] = A[i + sub_n][j]; A22[i][j] = A[i + sub_n][j + sub_n];
            B11[i][j] = B[i][j]; B12[i][j] = B[i][j + sub_n];
            B21[i][j] = B[i + sub_n][j]; B22[i][j] = B[i + sub_n][j + sub_n];
        }
    }
    
    // Calculate the 7 Products (P1 to P7)
    // P1 = A11 * (B12 - B22)
    subtract_matrices(sub_n, B12, B22, TempMatrix1);
    multiply_strassen_recursive(sub_n, A11, TempMatrix1, P1);
    
    // P2 = (A11 + A12) * B22
    add_matrices(sub_n, A11, A12, TempMatrix1);
    multiply_strassen_recursive(sub_n, TempMatrix1, B22, P2);
    
    // P3 = (A21 + A22) * B11
    add_matrices(sub_n, A21, A22, TempMatrix1);
    multiply_strassen_recursive(sub_n, TempMatrix1, B11, P3);
    
    // P4 = A22 * (B21 - B11)
    subtract_matrices(sub_n, B21, B11, TempMatrix1);
    multiply_strassen_recursive(sub_n, A22, TempMatrix1, P4);
    
    // P5 = (A11 + A22) * (B11 + B22)
    add_matrices(sub_n, A11, A22, TempMatrix1);
    add_matrices(sub_n, B11, B22, TempMatrix2);
    multiply_strassen_recursive(sub_n, TempMatrix1, TempMatrix2, P5);
    
    // P6 = (A12 - A22) * (B21 + B22)
    subtract_matrices(sub_n, A12, A22, TempMatrix1);
    add_matrices(sub_n, B21, B22, TempMatrix2);
    multiply_strassen_recursive(sub_n, TempMatrix1, TempMatrix2, P6);
    
    // P7 = (A11 - A21) * (B11 + B12)
    subtract_matrices(sub_n, A11, A21, TempMatrix1);
    add_matrices(sub_n, B11, B12, TempMatrix2);
    multiply_strassen_recursive(sub_n, TempMatrix1, TempMatrix2, P7);
    
    // Combine P's to get result sub-matrices C11, C12, C21, C22
    // C11 = P5 + P4 - P2 + P6
    add_matrices(sub_n, P5, P4, TempMatrix1);
    subtract_matrices(sub_n, TempMatrix1, P2, TempMatrix2);
    add_matrices(sub_n, TempMatrix2, P6, C11_res);
    
    // C12 = P1 + P2
    add_matrices(sub_n, P1, P2, C12_res);
    
    // C21 = P3 + P4
    add_matrices(sub_n, P3, P4, C21_res);
    
    // C22 = P5 + P1 - P3 - P7
    add_matrices(sub_n, P5, P1, TempMatrix1);
    subtract_matrices(sub_n, TempMatrix1, P3, TempMatrix2);
    subtract_matrices(sub_n, TempMatrix2, P7, C22_res);
    
    // Recombine C sub-matrices into final result C
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            C[i][j] = C11_res[i][j];
            C[i][j + sub_n] = C12_res[i][j];
            C[i + sub_n][j] = C21_res[i][j];
            C[i + sub_n][j + sub_n] = C22_res[i][j];
        }
    }
    
    // Clean up all 19 matrices
    destroy_matrix(sub_n, A11); destroy_matrix(sub_n, A12); destroy_matrix(sub_n, A21); destroy_matrix(sub_n, A22);
    destroy_matrix(sub_n, B11); destroy_matrix(sub_n, B12); destroy_matrix(sub_n, B21); destroy_matrix(sub_n, B22);
    destroy_matrix(sub_n, C11_res); destroy_matrix(sub_n, C12_res); destroy_matrix(sub_n, C21_res); destroy_matrix(sub_n, C22_res);
    destroy_matrix(sub_n, P1); destroy_matrix(sub_n, P2); destroy_matrix(sub_n, P3); destroy_matrix(sub_n, P4);
    destroy_matrix(sub_n, P5); destroy_matrix(sub_n, P6); destroy_matrix(sub_n, P7);
    destroy_matrix(sub_n, TempMatrix1); destroy_matrix(sub_n, TempMatrix2);
}

// --- Main Benchmark Function ---

int main() {
    int matrix_sizes[6] = {2, 4, 8, 16, 32, 64};
    const int benchmark_iterations = 1000;

    srand(time(NULL));
    
    printf("Strassen Matrix Multiplication Benchmark\n");
    printf("Matrix Size\tTotal Time (seconds)\n");
    printf("--------------------------------------\n");

    for(int i = 0; i < 6; i++){
        int current_size = matrix_sizes[i];
        
        double total_elapsed_time = 0.0;
        
        for(int iteration_count = 0; iteration_count < benchmark_iterations; iteration_count++) {
            int **MatrixA = create_square_matrix(current_size);
            int **MatrixB = create_square_matrix(current_size);
            int **MatrixC = create_square_matrix(current_size);
            
            for (int row = 0; row < current_size; row++) {
                for (int col = 0; col < current_size; col++) {
                    MatrixA[row][col] = rand() % 100;
                    MatrixB[row][col] = rand() % 100;
                }
            }
            
            clock_t start_time = clock();
            multiply_strassen_recursive(current_size, MatrixA, MatrixB, MatrixC);
            clock_t end_time = clock();
            
            total_elapsed_time += ((double)(end - start_time)) / CLOCKS_PER_SEC;
            
            destroy_matrix(current_size, MatrixA);
            destroy_matrix(current_size, MatrixB);
            destroy_matrix(current_size, MatrixC);
        }
        
        printf("%dx%d\t\t%lf\n", current_size, current_size, total_elapsed_time);
    }

    return 0;
}
