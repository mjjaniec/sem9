

void compute_parallel_simple(Star * stars, int N, int world_rank, int world_size, MPI_Datatype StarDataType) {
    int per_host = ceil (((double)N)/world_size);
    Star* local = (Star*) malloc(per_host * sizeof(Star));
    Star* remote = (Star*) malloc(per_host * sizeof(Star));

    MPI_Scatter(stars, per_host, StarDataType, local, per_host, StarDataType, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(stars, per_host, StarDataType, remote, per_host, StarDataType, ROOT, MPI_COMM_WORLD);

    // last node would probably get incomplete table, last elements should be nulled

    if (world_rank == world_size -1) {
        for (int i = N; i < per_host * world_size; ++i) {
            local[i].mass = remote[i].mass = 0.0;
        }
    }

    for (int iteration = 0;; ++iteration) {
        for (int i = 0; i < per_host; ++i) {
            for (int j = 0; j <per_host; ++j) {
                V3 diff = V3_sub(remote[j].position, local[i].position);
                if (V3_isZero(diff)) {
                    continue;
                }
                double coefficient = GRAVITY / math_qube(V3_norm(diff));
                local[i].acceleration = V3_add(local[i].acceleration, V3_mul_sv(coefficient * remote[j].mass, diff));
                remote[j].acceleration = V3_sub(remote[j].acceleration, V3_mul_sv(coefficient * local[i].mass, diff));
            }
        }
        if (iteration == world_size - 1) {
            break;
        }

        MPI_Send(remote, per_host, StarDataType, next(world_rank, world_size), TAG, MPI_COMM_WORLD);
        MPI_Recv(remote, per_host, StarDataType, prev(world_rank, world_size), TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Gather(local, per_host, StarDataType, stars, per_host, StarDataType, ROOT, MPI_COMM_WORLD);
}
