#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- Matrix Memory Management ---

int** create_square_matrix(int dimension_n) {
    int** matrix = (int**)malloc(dimension_n * sizeof(int*));
    for (int i = 0; i < dimension_n; i++) {
        matrix[i] = (int*)malloc(dimension_n * sizeof(int));
    }
    return matrix;
}

void destroy_matrix(int dimension_n, int** matrix) {
    for (int i = 0; i < dimension_n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// --- Standard O(n^3) Multiplication ---

void multiply_matrices_standard(int n, int** MatrixA, int** MatrixB, int** MatrixC) {
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            MatrixC[i][j] = 0;
            for(int k = 0; k < n; k++){
                MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
}

// --- Main Benchmark Program ---

int main() {
    int matrix_sizes[6] = {2, 4, 8, 16, 32, 64};
    const int benchmark_iterations = 1000;
    
    srand(time(NULL));

    printf("Matrix Size\tTotal Time (seconds)\n");
    printf("--------------------------------------\n");
    
    for(int i = 0; i < 6; i++){
        int current_size = matrix_sizes[i];
        
        double total_elapsed_time = 0.0;
        
        for(int iteration_count = 0; iteration_count < benchmark_iterations; iteration_count++) {
            int **MatrixA = create_square_matrix(current_size);
            int **MatrixB = create_square_matrix(current_size);
            int **MatrixC = create_square_matrix(current_size);
            
            // Populate matrices with random data
            for (int row = 0; row < current_size; row++) {
                for (int col = 0; col < current_size; col++) {
                    MatrixA[row][col] = rand() % 1000;
                    MatrixB[row][col] = rand() % 1000;
                }
            }
            
            clock_t start_time = clock();
            multiply_matrices_standard(current_size, MatrixA, MatrixB, MatrixC);
            clock_t end_time = clock();
            
            total_elapsed_time += ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
            
            destroy_matrix(current_size, MatrixA);
            destroy_matrix(current_size, MatrixB);
            destroy_matrix(current_size, MatrixC);
        }
        
        printf("%dx%d\t\t%lf\n", current_size, current_size, total_elapsed_time);
    }
    return 0;
}
