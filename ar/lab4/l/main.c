#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>
#include <stdbool.h>

const int PREFIX = 3;
const bool CIRCLE = false;
#include "graph.c"
const int ROOT = 0;


void brodcast_graph(Graph* graph, int world_rank) {
    void * bytes = world_rank == ROOT ? Graph_toBytes(graph) : malloc(Graph_bytesSize(graph));
    MPI_Bcast(bytes, Graph_bytesSize(graph), MPI_BYTE, ROOT, MPI_COMM_WORLD);
    if (world_rank != ROOT) {
        Graph_fromBytes(graph, bytes);
    }
    free(bytes);
    MPI_Barrier(MPI_COMM_WORLD);
}


double recurse(Graph* graph, bool* available, double actual, double best, int depth, int* order, int * bestOrder) {
    if (actual + graph->potentials[graph->size - depth + 1] > best) {
        return best;
    }
    if (depth == graph->size) {
        actual += Graph_cost(graph, order[0], order[depth-1]);
        if (actual < best) {
            memcpy(bestOrder, order, graph->size * sizeof(int));
            best = actual;
        }
        return best;
    } else {
        double result;
        for(int i = 0; i < graph->size; ++i) {
            if (available[i]) {
                double cost = actual + Graph_cost(graph, i, order[depth-1]);
                available[i] = false;
                order[depth] = i;
                result = recurse(graph, available, cost, best, depth + 1, order, bestOrder);
                if (result < best) {
                    best = result;
                }
                available[i] = true;
            }
        }
        return best;
    }
}


void master(Graph* graph, SalesmanPath* bestPath, int prefix, int world_size) {
    MPI_Status status;
    int all = perm_compute_all(graph->size - 1, prefix - 1);

    int data_size = sizeof(double) + sizeof(int) * graph->size;
    void* data = malloc(data_size);
    int* order = (int*)(data + sizeof(double));
    double* best = (double*)data;
   // bestPath->cost = 1000;

    for (int i = 0; i < all; ++i) {
        MPI_Recv(data, data_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // process results
        if (*best < bestPath->cost) {
            bestPath->cost = *best;
            memcpy(bestPath->vertexes, order, sizeof(int) * graph->size);
            printf("master best update: %f\n", (float)*best);
        }

        //create new task
        identity(order, graph->size);
        set_perm(order, graph->size, i, prefix);
        *best = bestPath->cost;

        MPI_Send(data, data_size, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
    }

    //stop workers
    for (int i = 1; i < world_size; ++i) {
        MPI_Recv(data, data_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (*best < bestPath->cost) {
            bestPath->cost = *best;
            memcpy(bestPath->vertexes, order, sizeof(int) * graph->size);
        }
        *best = 0;
        MPI_Send(data, data_size, MPI_BYTE, i, 0, MPI_COMM_WORLD);
    }

    free(data);
}


void worker(Graph* graph, int prefix) {
    int data_size = sizeof(double) + sizeof(int) * graph->size;
    void* data = malloc(data_size);
    double* best = ((double*)data);
    int* order = (int*)(data + sizeof(double));
    int* bestOrder = malloc(graph->size * sizeof(int));
    bool* available = malloc(graph->size * sizeof(bool));
    double actual;

    *best = 1.e10;      //value big enough to never be stored

    while (true) {
        //request for task
        MPI_Send(data, data_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(data, data_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        printf("%f\n", (float)(*best));
//        for (int i = 0; i < graph->size; ++i) {
//            printf("%d ", order[i]);
//        }
//        printf("\n");
        if (*best == 0) {       //stop condition
     //      printf("stopping worker\n");
            free(available);
            free(bestOrder);
            return;
        }
        // computation
        actual = 0;
        memset(available, true, graph->size);
        for (int i = 0; i < prefix; ++i) {
            available[order[i]] = false;
            if ( i > 0) {
                actual += Graph_cost(graph, order[i-1], order[i]);
            }
        }

        *best = recurse(graph, available, actual, *best, prefix, order, bestOrder);
        memcpy(order, bestOrder, graph->size * sizeof(int));
    }
}


void worker_sequent(Graph* graph, void* data, int prefix) {
    double* best = ((double*)data);
    int* order = (int*)(data + sizeof(double));
    int* bestOrder = malloc(graph->size * sizeof(int));
    bool* available = malloc(graph->size * sizeof(bool));
    double actual = 0;
    memset(available, true, graph->size);

    for (int i = 0; i < prefix; ++i) {
        available[order[i]] = false;
        if ( i > 0) {
            actual += Graph_cost(graph, order[i-1], order[i]);
        }
    }
    *best= recurse(graph, available, actual, *best, prefix, order, bestOrder);
    memcpy(data + sizeof(double), bestOrder, graph->size*sizeof(int));
    free(available);
    free(bestOrder);
}


void sequent(Graph* graph, SalesmanPath* bestPath, int prefix) {
    int all = perm_compute_all(graph->size - 1, prefix - 1);

    int data_size = sizeof(double) + sizeof(int) * graph->size;
    void* data = malloc(data_size);
    int* order = (int*)(data + sizeof(double));
    double* best = (double*)data;

    for (int i = 0; i < all; ++i) {
        identity(order, graph->size);
        set_perm(order, graph->size, i, prefix);
        *best = bestPath->cost;


        worker_sequent(graph, data, prefix);
        double cost = *best;
        if (cost < bestPath->cost) {
            bestPath->cost = cost;
            memcpy(bestPath->vertexes, data + sizeof(double), sizeof(int) * graph->size);
        }
    }

    free(data);
}


void compute_salesman(Graph* graph, SalesmanPath* bestPath, int world_rank, int world_size) {
    int prefix = PREFIX;
    brodcast_graph(graph, world_rank);
    bool seq = 0;
    if (seq) {
        sequent(graph, bestPath, prefix);
    } else {
        if (world_rank == ROOT) {
            master(graph, bestPath, prefix, world_size);
        } else {
            worker(graph, prefix);
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
           if (world_rank == 0) {
            printf("Usage: %s #cities\n", argv[0]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return -1;
    }


    int N = atoi(argv[1]);
    srand(0);

    Graph graph;
    SalesmanPath bestPath;
    if (world_rank == ROOT) {
        Graph_init(&graph, N);
        bestPath = Graph_getSimpleSalesmanPath(&graph);
    } else {
        // just to properly compute graph size
        graph.size = N;
    }

    struct timeval  tv0, tv1;
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == ROOT) {
        gettimeofday(&tv0, NULL);
    }

    compute_salesman(&graph, &bestPath, world_rank, world_size);

    if (world_rank == ROOT) {
        gettimeofday(&tv1, NULL);
        float time = tv1.tv_sec - tv0.tv_sec + (tv1.tv_usec - tv0.tv_usec) / 1000000.0f;

        printf("N: %d, time: %f\n", N, time);
        printf("Cost: %f, order: ", bestPath.cost);
        for (int i = 0; i <bestPath.size; ++i) {
            printf("%d ", bestPath.vertexes[i]);
        }
        printf("\n");
    }
    MPI_Finalize();
    return 0;
}
