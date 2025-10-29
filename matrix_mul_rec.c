#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- Matrix Memory and Utility ---

int** create_square_matrix(int dimension_n) {
    int** matrix_handle = (int**)malloc(dimension_n * sizeof(int*));
    for (int i = 0; i < dimension_n; i++) {
        matrix_handle[i] = (int*)malloc(dimension_n * sizeof(int));
    }
    return matrix_handle;
}

void destroy_matrix(int** matrix_handle, int dimension_n) {
    for (int i = 0; i < dimension_n; i++) {
        free(matrix_handle[i]);
    }
    free(matrix_handle);
}

void add_matrices(int** MatrixA, int** MatrixB, int** MatrixResult, int dimension_n) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            MatrixResult[i][j] = MatrixA[i][j] + MatrixB[i][j];
        }
    }
}

void display_matrix(int** matrix, int dimension_n) {
    for (int i = 0; i < dimension_n; i++) {
        for (int j = 0; j < dimension_n; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// --- Divide and Conquer O(n^3) Multiplication ---

void multiply_matrices_recursively(int** A, int** B, int** C, int n) {
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    
    if (n == 2) {
        C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
        C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
        C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
        C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
        return;
    }
    
    int sub_dimension = n / 2;
    
    // Allocate 13 temporary matrices for sub-problems
    int** A_top_left = create_square_matrix(sub_dimension);
    int** A_top_right = create_square_matrix(sub_dimension);
    int** A_bottom_left = create_square_matrix(sub_dimension);
    int** A_bottom_right = create_square_matrix(sub_dimension);
    
    int** B_top_left = create_square_matrix(sub_dimension);
    int** B_top_right = create_square_matrix(sub_dimension);
    int** B_bottom_left = create_square_matrix(sub_dimension);
    int** B_bottom_right = create_square_matrix(sub_dimension);
    
    int** C_top_left = create_square_matrix(sub_dimension);
    int** C_top_right = create_square_matrix(sub_dimension);
    int** C_bottom_left = create_square_matrix(sub_dimension);
    int** C_bottom_right = create_square_matrix(sub_dimension);
    
    int** TempStorage = create_square_matrix(sub_dimension);

    // Partition A and B
    for (int i = 0; i < sub_dimension; i++) {
        for (int j = 0; j < sub_dimension; j++) {
            A_top_left[i][j] = A[i][j];
            A_top_right[i][j] = A[i][j + sub_dimension];
            A_bottom_left[i][j] = A[i + sub_dimension][j];
            A_bottom_right[i][j] = A[i + sub_dimension][j + sub_dimension];
            
            B_top_left[i][j] = B[i][j];
            B_top_right[i][j] = B[i][j + sub_dimension];
            B_bottom_left[i][j] = B[i + sub_dimension][j];
            B_bottom_right[i][j] = B[i + sub_dimension][j + sub_dimension];
        }
    }
    
    // Calculate C11 = A11*B11 + A12*B21
    multiply_matrices_recursively(A_top_left, B_top_left, C_top_left, sub_dimension);
    multiply_matrices_recursively(A_top_right, B_bottom_left, TempStorage, sub_dimension);
    add_matrices(C_top_left, TempStorage, C_top_left, sub_dimension);
    
    // Calculate C12 = A11*B12 + A12*B22
    multiply_matrices_recursively(A_top_left, B_top_right, C_top_right, sub_dimension);
    multiply_matrices_recursively(A_top_right, B_bottom_right, TempStorage, sub_dimension);
    add_matrices(C_top_right, TempStorage, C_top_right, sub_dimension);
    
    // Calculate C21 = A21*B11 + A22*B21
    multiply_matrices_recursively(A_bottom_left, B_top_left, C_bottom_left, sub_dimension);
    multiply_matrices_recursively(A_bottom_right, B_bottom_left, TempStorage, sub_dimension);
    add_matrices(C_bottom_left, TempStorage, C_bottom_left, sub_dimension);
    
    // Calculate C22 = A21*B12 + A22*B22
    multiply_matrices_recursively(A_bottom_left, B_top_right, C_bottom_right, sub_dimension);
    multiply_matrices_recursively(A_bottom_right, B_bottom_right, TempStorage, sub_dimension);
    add_matrices(C_bottom_right, TempStorage, C_bottom_right, sub_dimension);
    
    // Combine sub-results into final matrix C
    for (int i = 0; i < sub_dimension; i++) {
        for (int j = 0; j < sub_dimension; j++) {
            C[i][j] = C_top_left[i][j];
            C[i][j + sub_dimension] = C_top_right[i][j];
            C[i + sub_dimension][j] = C_bottom_left[i][j];
            C[i + sub_dimension][j + sub_dimension] = C_bottom_right[i][j];
        }
    }
    
    // Clean up temporary memory
    destroy_matrix(A_top_left, sub_dimension); destroy_matrix(A_top_right, sub_dimension);
    destroy_matrix(A_bottom_left, sub_dimension); destroy_matrix(A_bottom_right, sub_dimension);
    destroy_matrix(B_top_left, sub_dimension); destroy_matrix(B_top_right, sub_dimension);
    destroy_matrix(B_bottom_left, sub_dimension); destroy_matrix(B_bottom_right, sub_dimension);
    destroy_matrix(C_top_left, sub_dimension); destroy_matrix(C_top_right, sub_dimension);
    destroy_matrix(C_bottom_left, sub_dimension); destroy_matrix(C_bottom_right, sub_dimension);
    destroy_matrix(TempStorage, sub_dimension); 
}

// --- Main Benchmark Function ---

int main() {
    int matrix_sizes[6] = {2, 4, 8, 16, 32, 64};
    const int benchmark_iterations = 1000;
    
    srand(time(NULL));
    
    printf("Recursive Matrix Multiplication Benchmark\n");
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
                    MatrixC[row][col] = 0;
                }
            }
            
            clock_t start_time = clock();
            multiply_matrices_recursively(MatrixA, MatrixB, MatrixC, current_size);
            clock_t end_time = clock();
            
            total_elapsed_time += ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
            
            destroy_matrix(MatrixA, current_size);
            destroy_matrix(MatrixB, current_size);
            destroy_matrix(MatrixC, current_size);
        }
        
        printf("%dx%d\t\t%lf\n", current_size, current_size, total_elapsed_time);
    }
    return 0;
}
