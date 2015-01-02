#include <mpi/mpi.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#define TAG 0
#define ROOT 0
#define $ printf("Line: %d\n",__LINE__);


void init_master_matrix(int N, double **M) {
    *M = (double*) malloc(N * (N + 15) * sizeof(double));

    double x;
    for (int row = 0; row < N; ++row) {
        for (int column = 0; column < N; ++column) {
            x = row / (N - 1.0);
            (*M)[row * N + column] = 1.0 - 4 * (x - 0.5) *  (x - 0.5);
        }
    }
}

void init_worker_matrix(int N, double **M, int world_size) {
    int part = (N -2) / world_size + 3;

    *M = (double*) malloc(N * part * sizeof(double));
}

void free_matrix(double *M) {
    free(M);
}

int bottom_edge(int N, int world_size) {
    int part = (N - 2) / world_size;
    int ret = 0;
    if ((N - 2) %world_size != 0) {
        part += 1;
        ret = part + (N - (world_size * part + 2));
        ret += 1;
    }
    else {
        ret = part +1;
    }
    return ret;
}

void print_matrix(int world_rank, double * matrix, int rows, int columns) {
    printf("==== %d ====\n", world_rank);
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            printf("%6.2f", (float) matrix[row * columns + column]);
        }
        printf("\n");
    }
    printf("\n");
}

void do_compute(double *W, int N, int part, int world_rank, int world_size) {
//    sleep(1*world_rank);
//    print_matrix(world_rank, W, part+2, N);
    int init_row = 1;

    int max_row = part + 1;
    if (world_rank == world_size - 1) {
        max_row = bottom_edge(N, world_size);
    }
    for (int row = init_row; row < max_row; ++row) {
        for (int column = 1; column < N - 1; ++column) {
            int index = N * row + column;
            W[index] = 0.25 * (W[index-N] + W[index-1] + W[index+1] + W[index+N]);
        }
    }
//    print_matrix(world_rank, W, part+2, N);

}

void compute(int world_rank, int world_size, int iterations) {
    int k;
    for(int iter = 0; iter<iterations; ++iter) {
        if (world_rank == 0) {
            MPI_Send(&k, 1, MPI_BYTE,1,TAG,MPI_COMM_WORLD);
            MPI_Recv(&k, 1, MPI_BYTE,1,TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&k, 1, MPI_BYTE,0,TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&k, 1, MPI_BYTE,0,TAG,MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char** argv) {
    int world_rank;
    int world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 2) {
        printf("usage: %s problem_size\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    srand(time(NULL));
    double *M;
    double *W;
    int N = atoi(argv[1]);

    bottom_edge(N, world_size);

    if (world_rank == 0) {
        init_master_matrix(N, &M);
    }

    init_worker_matrix(N, &W, world_size);

   // print_matrix(M, N, N);

    struct timeval  tv0, tv1;
    if (world_rank == 0) {
        gettimeofday(&tv0, NULL);
    }

    compute(world_rank, world_size, 1);

    if (world_rank == 0) {
        gettimeofday(&tv1, NULL);
//        print_matrix(0, M, N, N);
        float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;

        printf("N: %d, time: %f\n", N, time);
        free_matrix(M);
    }
    free_matrix(W);
    MPI_Finalize();
    return 0;
}
