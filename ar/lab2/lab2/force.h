#define TAG 0
#define ROOT 0
#define $ if(world_rank == 0) {printf("Line: %d\n",__LINE__);}

int next(int world_rank, int world_size, int count) {
    return (world_rank + count) % world_size;
}

int prev(int world_rank, int world_size, int count) {
    return (world_rank + world_size - count) % world_size;
}

void compute_force_sequential(Star * stars, int N) {
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

void compute_force_parallel(Star * stars, int N, int world_rank, int world_size, MPI_Datatype StarDataType) {
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
