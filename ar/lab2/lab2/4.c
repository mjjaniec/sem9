#include <mpi/mpi.h>

#include "v3.h"
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#define TAG 0
#define ROOT 0
#define $ if(world_rank == 0) {printf("Line: %d\n",__LINE__);}

int next(int world_rank, int world_size, int count) {
    return (world_rank + count) % world_size;
}

int prev(int world_rank, int world_size, int count) {
    return (world_rank + world_size - count) % world_size;
}

void compute_sequential(Star * stars, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j){
                V3 diff = V3_sub(stars[j].position, stars[i].position);
                double scalar = GRAVITY * stars[i].mass * stars[j].mass / math_qube(V3_norm(diff));
                stars[i].force = V3_add(stars[i].force, V3_mul_sv(scalar, diff));
            }
        }
    }
}

void compute_part(Star* local, Star* remote, int per_host, bool accumulate) {
    for (int i = 0; i < per_host; ++i) {
        for (int j = 0; j <per_host; ++j) {
            V3 diff = V3_sub(remote[j].position, local[i].position);
            if (V3_isZero(diff)) {
                continue;
            }
            double scalar = GRAVITY * local[i].mass * remote[j].mass / math_qube(V3_norm(diff));
            local[i].force = V3_add(local[i].force, V3_mul_sv(scalar, diff));

            if (accumulate) {
                //accumulator
                remote[j].force = V3_sub(remote[j].force, V3_mul_sv(scalar, diff));
            }
        }
    }
}

void use_accumulators(Star* local, Star* remote, Star* remote_rec,
                      int per_host, int world_rank, int world_size, int step, MPI_Datatype StarDataType) {
    MPI_Sendrecv(remote,     per_host, StarDataType, prev(world_rank, world_size, step), TAG,
                 remote_rec, per_host, StarDataType, next(world_rank, world_size, step), TAG,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(remote, remote_rec, per_host * sizeof(Star));
    for (int i = 0; i<per_host; ++i) {
        local[i].force = V3_add(local[i].force, remote[i].force);
    }
}

void compute_parallel(Star * stars, int N, int world_rank, int world_size, MPI_Datatype StarDataType) {
    //number of stars per node
    int per_host = ceil (((double)N)/world_size);
    //stars kept on node
    Star* local = (Star*) malloc(per_host * sizeof(Star));
    //remote stars - thier force field is unused, would be used as accumulator
    Star* remote = (Star*) malloc(per_host * sizeof(Star));
    Star* remote_rec = (Star*) malloc(per_host * sizeof(Star));

    MPI_Scatter(stars, per_host, StarDataType, local, per_host, StarDataType, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(stars, per_host, StarDataType, remote, per_host, StarDataType, ROOT, MPI_COMM_WORLD);

    // last node would probably get incomplete table, last elements mass set to null would eliminate thier impact
    if (per_host * world_size > N && world_rank == world_size -1) {
        for (int i = N - per_host * (world_size - 1); i < per_host; ++i) {
            local[i].mass = remote[i].mass = 0.0;
        }
    }

    int last_step = (world_size + 1) / 2;

    for (int step = 0; step < last_step; ++step) {
        compute_part(local, remote, per_host, step > 0);
        MPI_Sendrecv(remote,     per_host, StarDataType, next(world_rank, world_size, 1), TAG,
                     remote_rec, per_host, StarDataType, prev(world_rank, world_size, 1), TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        memcpy(remote, remote_rec, per_host * sizeof(Star));
    }

    if (world_size % 2 == 0) {
        compute_part(local, remote, per_host, false);
    }
    if (world_size > 2) {
        use_accumulators(local, remote, remote_rec, per_host, world_rank, world_size, last_step, StarDataType);\
    }

    MPI_Gather(local, per_host, StarDataType, stars, per_host, StarDataType, ROOT, MPI_COMM_WORLD);
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
//    V3*sequence = NULL;

    if (world_rank == ROOT) {
        init_stars(&stars, N, world_size);
//        sequence = (V3*)malloc(N * sizeof(V3));

    }

//    MPI_Barrier(MPI_COMM_WORLD);
//    if (world_rank == ROOT) {
//        printf("sequential: \n");

//        gettimeofday(&tv0, NULL);

//        compute_sequential(stars, N);

//        gettimeofday(&tv1, NULL);
//        float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;

//        printf("N: %d, time: %f\n", N, time);

//        for (int i = 0; i<N; ++i) {
//            sequence[i] = stars[i].force;
//            stars[i].force = V3_create(0,0,0);
////            V3_print(sequence[i], true);
//        }

//        printf("\n===============\n\nparallel:\n");
//    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (world_rank == ROOT) {
        gettimeofday(&tv0, NULL);
    }


    compute_parallel(stars, N, world_rank, world_size, StarDataType);


    if (world_rank == ROOT) {
        gettimeofday(&tv1, NULL);
        float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;

        printf("N: %d, time: %f\n", N, time);
//        for (int i = 0; i<N; ++i) {
//            if (!V3_equals(sequence[i], stars[i].force)) {
//                printf("%d:\t", i);
//                V3_print(V3_sub(sequence[i], stars[i].force), true);
//            }
//        }
    }
    MPI_Finalize();
    return 0;
}
