#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- Core Matrix Utilities ---

int** allocate_square_matrix(int dimension_n) {
    int** matrix = (int**)malloc(dimension_n * sizeof(int*));
    for (int i = 0; i < dimension_n; i++) {
        matrix[i] = (int*)malloc(dimension_n * sizeof(int));
    }
    return matrix;
}

void release_matrix_memory(int dimension_n, int** matrix) {
    for (int i = 0; i < dimension_n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void add_matrices(int** MatrixA, int** MatrixB, int** MatrixResult, int dimension_n) {
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

void print_matrix_content(int** matrix, int dimension_n) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// --- Standard O(n^3) Algorithm ---

void multiply_standard(int n, int** MatrixA, int** MatrixB, int** MatrixC) {
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            MatrixC[i][j] = 0;
            for(int k = 0; k < n; k++){
                MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
}

// --- Simple Divide and Conquer O(n^3) Algorithm ---

void multiply_divide_and_conquer(int** A, int** B, int** C, int n) {
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    
    int sub_n = n / 2;
    
    // Allocate 10 temporary matrices (sub-matrices of A, B, C, and one temporary for addition)
    int** A11 = allocate_square_matrix(sub_n);
    int** A12 = allocate_square_matrix(sub_n);
    int** A21 = allocate_square_matrix(sub_n);
    int** A22 = allocate_square_matrix(sub_n);
    
    int** B11 = allocate_square_matrix(sub_n);
    int** B12 = allocate_square_matrix(sub_n);
    int** B21 = allocate_square_matrix(sub_n);
    int** B22 = allocate_square_matrix(sub_n);
    
    int** C11_temp = allocate_square_matrix(sub_n); // C11 = A11*B11 + A12*B21
    int** C12_temp = allocate_square_matrix(sub_n);
    int** C21_temp = allocate_square_matrix(sub_n);
    int** C22_temp = allocate_square_matrix(sub_n);

    int** TempStorage = allocate_square_matrix(sub_n); // Used for A12*B21, etc.
    
    // Partition A and B into sub-matrices
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + sub_n];
            A21[i][j] = A[i + sub_n][j];
            A22[i][j] = A[i + sub_n][j + sub_n];
            
            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + sub_n];
            B21[i][j] = B[i + sub_n][j];
            B22[i][j] = B[i + sub_n][j + sub_n];
        }
    }
    
    // C11 = A11*B11 + A12*B21
    multiply_divide_and_conquer(A11, B11, C11_temp, sub_n); // A11*B11 stored in C11_temp
    multiply_divide_and_conquer(A12, B21, TempStorage, sub_n); // A12*B21 stored in TempStorage
    add_matrices(C11_temp, TempStorage, C11_temp, sub_n); // C11_temp = C11_temp + TempStorage
    
    // C12 = A11*B12 + A12*B22
    multiply_divide_and_conquer(A11, B12, C12_temp, sub_n);
    multiply_divide_and_conquer(A12, B22, TempStorage, sub_n);
    add_matrices(C12_temp, TempStorage, C12_temp, sub_n);
    
    // C21 = A21*B11 + A22*B21
    multiply_divide_and_conquer(A21, B11, C21_temp, sub_n);
    multiply_divide_and_conquer(A22, B21, TempStorage, sub_n);
    add_matrices(C21_temp, TempStorage, C21_temp, sub_n);
    
    // C22 = A21*B12 + A22*B22
    multiply_divide_and_conquer(A21, B12, C22_temp, sub_n);
    multiply_divide_and_conquer(A22, B22, TempStorage, sub_n);
    add_matrices(C22_temp, TempStorage, C22_temp, sub_n);
    
    // Combine result sub-matrices into C
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            C[i][j] = C11_temp[i][j];
            C[i][j + sub_n] = C12_temp[i][j];
            C[i + sub_n][j] = C21_temp[i][j];
            C[i + sub_n][j + sub_n] = C22_temp[i][j];
        }
    }
    
    // Release all temporary memory
    release_matrix_memory(sub_n, A11); release_matrix_memory(sub_n, A12);
    release_matrix_memory(sub_n, A21); release_matrix_memory(sub_n, A22);
    release_matrix_memory(sub_n, B11); release_matrix_memory(sub_n, B12);
    release_matrix_memory(sub_n, B21); release_matrix_memory(sub_n, B22);
    release_matrix_memory(sub_n, C11_temp); release_matrix_memory(sub_n, C12_temp);
    release_matrix_memory(sub_n, C21_temp); release_matrix_memory(sub_n, C22_temp);
    release_matrix_memory(sub_n, TempStorage);
}

// --- Strassen's O(n^2.807) Algorithm ---

void multiply_strassen(int n, int **A, int **B, int **C) {
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    
    int sub_n = n / 2;
    
    // Allocate 19 matrices: 4 for A, 4 for B, 4 for C, 7 for P, and 2 temps
    int **A11 = allocate_square_matrix(sub_n); int **A12 = allocate_square_matrix(sub_n);
    int **A21 = allocate_square_matrix(sub_n); int **A22 = allocate_square_matrix(sub_n);
    int **B11 = allocate_square_matrix(sub_n); int **B12 = allocate_square_matrix(sub_n);
    int **B21 = allocate_square_matrix(sub_n); int **B22 = allocate_square_matrix(sub_n);
    int **C11 = allocate_square_matrix(sub_n); int **C12 = allocate_square_matrix(sub_n);
    int **C21 = allocate_square_matrix(sub_n); int **C22 = allocate_square_matrix(sub_n);
    int **P1 = allocate_square_matrix(sub_n); int **P2 = allocate_square_matrix(sub_n);
    int **P3 = allocate_square_matrix(sub_n); int **P4 = allocate_square_matrix(sub_n);
    int **P5 = allocate_square_matrix(sub_n); int **P6 = allocate_square_matrix(sub_n);
    int **P7 = allocate_square_matrix(sub_n);
    int **TempAddition = allocate_square_matrix(sub_n); 
    int **TempSubtraction = allocate_square_matrix(sub_n);
    
    // Partition A and B
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            A11[i][j] = A[i][j]; A12[i][j] = A[i][j + sub_n];
            A21[i][j] = A[i + sub_n][j]; A22[i][j] = A[i + sub_n][j + sub_n];
            B11[i][j] = B[i][j]; B12[i][j] = B[i][j + sub_n];
            B21[i][j] = B[i + sub_n][j]; B22[i][j] = B[i + sub_n][j + sub_n];
        }
    }
    
    // Calculate the 7 products (P1 to P7)
    // P1 = A11 * (B12 - B22)
    subtract_matrices(sub_n, B12, B22, TempSubtraction);
    multiply_strassen(sub_n, A11, TempSubtraction, P1);
    
    // P2 = (A11 + A12) * B22
    add_matrices(A11, A12, TempAddition, sub_n);
    multiply_strassen(sub_n, TempAddition, B22, P2);
    
    // P3 = (A21 + A22) * B11
    add_matrices(A21, A22, TempAddition, sub_n);
    multiply_strassen(sub_n, TempAddition, B11, P3);
    
    // P4 = A22 * (B21 - B11)
    subtract_matrices(sub_n, B21, B11, TempSubtraction);
    multiply_strassen(sub_n, A22, TempSubtraction, P4);
    
    // P5 = (A11 + A22) * (B11 + B22)
    add_matrices(A11, A22, TempAddition, sub_n);
    add_matrices(B11, B22, TempSubtraction, sub_n);
    multiply_strassen(sub_n, TempAddition, TempSubtraction, P5);
    
    // P6 = (A12 - A22) * (B21 + B22)
    subtract_matrices(sub_n, A12, A22, TempSubtraction);
    add_matrices(B21, B22, TempAddition, sub_n);
    multiply_strassen(sub_n, TempSubtraction, TempAddition, P6);
    
    // P7 = (A11 - A21) * (B11 + B12)
    subtract_matrices(sub_n, A11, A21, TempSubtraction);
    add_matrices(B11, B12, TempAddition, sub_n);
    multiply_strassen(sub_n, TempSubtraction, TempAddition, P7);
    
    // Combine P's to get result sub-matrices C11, C12, C21, C22
    // C11 = P5 + P4 - P2 + P6
    add_matrices(P5, P4, TempAddition, sub_n); // TempAddition = P5 + P4
    subtract_matrices(sub_n, TempAddition, P2, TempSubtraction); // TempSubtraction = P5 + P4 - P2
    add_matrices(TempSubtraction, P6, C11, sub_n); // C11 = TempSubtraction + P6
    
    // C12 = P1 + P2
    add_matrices(P1, P2, C12, sub_n);
    
    // C21 = P3 + P4
    add_matrices(P3, P4, C21, sub_n);
    
    // C22 = P5 + P1 - P3 - P7
    add_matrices(P5, P1, TempAddition, sub_n); // TempAddition = P5 + P1
    subtract_matrices(sub_n, TempAddition, P3, TempSubtraction); // TempSubtraction = P5 + P1 - P3
    subtract_matrices(sub_n, TempSubtraction, P7, C22); // C22 = TempSubtraction - P7
    
    // Recombine C sub-matrices into final result C
    for (int i = 0; i < sub_n; i++) {
        for (int j = 0; j < sub_n; j++) {
            C[i][j] = C11[i][j];
            C[i][j + sub_n] = C12[i][j];
            C[i + sub_n][j] = C21[i][j];
            C[i + sub_n][j + sub_n] = C22[i][j];
        }
    }
    
    // Release all 19 temporary memory allocations
    release_matrix_memory(sub_n, A11); release_matrix_memory(sub_n, A12); release_matrix_memory(sub_n, A21); release_matrix_memory(sub_n, A22);
    release_matrix_memory(sub_n, B11); release_matrix_memory(sub_n, B12); release_matrix_memory(sub_n, B21); release_matrix_memory(sub_n, B22);
    release_matrix_memory(sub_n, C11); release_matrix_memory(sub_n, C12); release_matrix_memory(sub_n, C21); release_matrix_memory(sub_n, C22);
    release_matrix_memory(sub_n, P1); release_matrix_memory(sub_n, P2); release_matrix_memory(sub_n, P3); release_matrix_memory(sub_n, P4);
    release_matrix_memory(sub_n, P5); release_matrix_memory(sub_n, P6); release_matrix_memory(sub_n, P7);
    release_matrix_memory(sub_n, TempAddition); release_matrix_memory(sub_n, TempSubtraction);
}

// --- Main Benchmark Function ---

int main(){
    int matrix_sizes[6] = {2, 4, 8, 16, 32, 64};
    const int num_iterations = 1000;
    
    // File setup
    FILE *fp_standard = fopen("standard_results.txt", "w");
    FILE *fp_divideconquer = fopen("divideconquer_results.txt", "w");
    FILE *fp_strassen = fopen("strassen_results.txt", "w");
    
    if (!fp_standard || !fp_divideconquer || !fp_strassen) {
        printf("Error: Could not open one or more result files.\n");
        return 1;
    }
    
    fprintf(fp_standard, "size,time\n");
    fprintf(fp_divideconquer, "size,time\n");
    fprintf(fp_strassen, "size,time\n");
    
    srand(time(NULL));
    
    printf("--- Matrix Multiplication Benchmark ---\n");
    
    for(int i = 0; i < 6; i++){
        int n = matrix_sizes[i];
        printf("Testing matrix size: %dx%d (Iterations: %d)\n", n, n, num_iterations);
        
        // --- 1. Standard Multiplication (O(n^3)) ---
        double total_time_standard = 0.0;
        
        for(int iter = 0; iter < num_iterations; iter++) {
            int **A = allocate_square_matrix(n);
            int **B = allocate_square_matrix(n);
            int **C = allocate_square_matrix(n);
            
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    A[row][col] = rand() % 100;
                    B[row][col] = rand() % 100;
                }
            }
            
            clock_t start = clock();
            multiply_standard(n, A, B, C);
            clock_t end = clock();
            
            total_time_standard += ((double)(end - start)) / CLOCKS_PER_SEC;
            
            release_matrix_memory(n, A);
            release_matrix_memory(n, B);
            release_matrix_memory(n, C);
        }
        
        printf("Standard O(n^3) \t- Total time: %lf seconds\n", total_time_standard);
        fprintf(fp_standard, "%d,%lf\n", n, total_time_standard);
        
        // --- 2. Divide and Conquer Multiplication (O(n^3)) ---
        double total_time_dc = 0.0;
        
        for(int iter = 0; iter < num_iterations; iter++) {
            int **A = allocate_square_matrix(n);
            int **B = allocate_square_matrix(n);
            int **C = allocate_square_matrix(n);
            
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    A[row][col] = rand() % 100;
                    B[row][col] = rand() % 100;
                }
            }
            
            clock_t start = clock();
            multiply_divide_and_conquer(A, B, C, n);
            clock_t end = clock();
            
            total_time_dc += ((double)(end - start)) / CLOCKS_PER_SEC;
            
            release_matrix_memory(n, A);
            release_matrix_memory(n, B);
            release_matrix_memory(n, C);
        }
        
        printf("Divide & Conquer O(n^3) - Total time: %lf seconds\n", total_time_dc);
        fprintf(fp_divideconquer, "%d,%lf\n", n, total_time_dc);
        
        // --- 3. Strassen Multiplication (O(n^2.807)) ---
        double total_time_strassen = 0.0;
        
        for(int iter = 0; iter < num_iterations; iter++) {
            int **A = allocate_square_matrix(n);
            int **B = allocate_square_matrix(n);
            int **C = allocate_square_matrix(n);
            
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    A[row][col] = rand() % 100;
                    B[row][col] = rand() % 100;
                }
            }
            
            clock_t start = clock();
            multiply_strassen(n, A, B, C);
            clock_t end = clock();
            
            total_time_strassen += ((double)(end - start)) / CLOCKS_PER_SEC;
            
            release_matrix_memory(n, A);
            release_matrix_memory(n, B);
            release_matrix_memory(n, C);
        }
        
        printf("Strassen O(n^2.807) \t- Total time: %lf seconds\n", total_time_strassen);
        fprintf(fp_strassen, "%d,%lf\n", n, total_time_strassen);
        printf("---\n");
    }
    
    fclose(fp_standard);
    fclose(fp_divideconquer);
    fclose(fp_strassen);
    printf("Benchmark complete. Results saved to files for plotting.\n");
    return 0;
}
