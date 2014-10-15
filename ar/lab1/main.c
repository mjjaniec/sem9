#include <mpi/mpi.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


void init_matrix(int N, double **M) {
    *M = (double*) malloc(N * N * sizeof(double));

    double x;
    for (int row = 0; row < N; ++row) {
        for (int column = 0; column < N; ++column) {
            x = row / (N - 1.0);
            *M[row * N + column] = 1.0 - 4 * (x - 0.5) *  (x - 0.5);
        }
    }
}

void free_matrix(double *M) {
    free(M);
}

void print_matrix(double * matrix, int rows, int columns) {
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            printf("%4d", (int) matrix[row * columns + column]);
        }
        printf("\n");
    }
    printf("\n");
}

void compute(double *M, int N, int proc_a, int proc_b) {

}

int main(int argc, char** argv) {
    int world_rank;
    int world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 4) {
        printf("usage: %s problem_size a b\n", argv[0]);
        printf("   where number_of_processors = a * b\n");
        MPI_Finalize();
        return -1;
    }

  /*  srand(time(NULL));
    double *M;
    int N = atoi(argv[1]);
    int A = atoi(argv[2]);
    int B = atoi(argv[3]);

    init_matrix(N, &M);
    //print_matrix("A", A, N, N);
    //print_matrix("B", B, N, 1);

    struct timeval  tv0, tv1;
    gettimeofday(&tv0, NULL);

    compute(M, N, A, B);

    gettimeofday(&tv1, NULL);
    //print_matrix("C", C, N, 1);
    float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;
    printf("N: %d, time: %f\n", N, time);

    free_matrix(M);
    MPI_Finalize();
    return 0;*/
}
