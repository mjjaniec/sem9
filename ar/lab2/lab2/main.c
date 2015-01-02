#include <mpi/mpi.h>

#include "v3.h"
#include "force.h"
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>


void compute_trajectory(Star* stars, int N, double dt, double max_t,
                        int world_rank, int world_size, MIP_Datatype StarDataType) {
    V3* velocity = (V3*)malloc(N * sizeof(V3));
    memset(velocity, 0, N * sizeof(V3));

    compute_force_parallel(stars, N, world_rank, world_size, StarDataType);
    for (double t = 0; t < max_t; t+=dt) {
        for (int i = 0; i < N; ++i) {
            V3_add_assign(velocity[i], V3_mul_sv(0.5 * dt / stars[i].mass, stars[i].force));
            V3_add_assign(stars[i].position, V3_mul_sv(dt, velocity[i]));
        }
        compute_force_parallel(stars, N, world_rank, world_szie, StarDataType);
        for (int i = 0; i < N; ++i) {
            V3_add_assign(velocity[i], V3_mul_sv(0.5 * dt / stars[i].mass, stars[i].force));
        }
    }

}


void init_stars(Star**stars, int N, int world_size) {
    int part = ceil (((double)N)/world_size);
    *stars = (Star*)malloc(part * world_size * sizeof(Star));
    for (int i = 0; i < N; ++i) {
        Star_init(&((*stars)[i]));
    }
}


int main(int argc, char** argv) {
    int world_rank;
    int world_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc != 2) {
           if (world_rank == 0) {
            printf("Usage: %s #starts\n", argv[0]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return -1;
    }


    srand(time(NULL));

    struct timeval  tv0, tv1;
    MPI_Datatype StarDataType;
    MPI_Type_contiguous(sizeof(Star), MPI_BYTE, &StarDataType);
    MPI_Type_commit(&StarDataType);

    int N = atoi(argv[1]);

    Star* stars;

    if (world_rank == ROOT) {
        init_stars(&stars, N, world_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == ROOT) {
        gettimeofday(&tv0, NULL);
    }

    compute_trajectory(stars, N, world_rank, world_size, StarDataType);

    if (world_rank == ROOT) {
        gettimeofday(&tv1, NULL);
        float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;

        printf("N: %d, time: %f\n", N, time);
    }
    MPI_Finalize();
    return 0;
}
