#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "util.c"
#include "disjont_set.c"
#include "list.c"
#include "stdbool.h"

typedef struct __Graph {
    int size;
    double* cost;
} Graph;

typedef struct __SalesmanPath {
    int size;
    int *vertexes;
    double cost;
} SalesmanPath;

typedef struct __Edge {
    int a, b;
    double cost;
} Edge;

Edge Edge_create(int a, int b, double cost) {
    Edge edge;
    edge.a = a;
    edge.b = b;
    edge.cost = cost;
    return edge;
}

int Edge_compare(const void* a, const void*b) {
    return sign(((Edge*)a)->cost - ((Edge*)b)->cost);
}

void Graph_init(Graph* self, int size) {
    self->cost = (double*)malloc(size * size * sizeof(double));
    self->size = size;

    double * x = malloc(size * sizeof(double));
    double * y = malloc(size * sizeof(double));

    int* indexes = malloc(size * sizeof(double));
    identity(indexes, size);
    random_perm(indexes, size);

    for (int i = 0; i < size; ++i) {
//        x[i] = rand() % 25;
//        y[i] = rand() % 25;
        x[i] = 25 * sin(indexes[i]*2.0*3.14159265/size);
        y[i] = 25 * cos(indexes[i]*2.0*3.14159265/size);
    }
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dist = sqrt(dx * dx + dy * dy);
            self->cost[i * size + j] = dist;
        }
    }
    free(x);
    free(y);
}

void Graph_free(Graph* self) {
    free(self->cost);
}

double Graph_cost(Graph* self, int i, int j) {
    return self->cost[self->size * i + j] ;
}


int Graph_bytesSize(Graph* self) {
    return sizeof(int) + self->size * self->size * sizeof(double);
}

/**
 * @brief Graph_toBytes those bytes needs to be free later
 */
void* Graph_toBytes(Graph* self) {
    void* data = malloc(Graph_bytesSize(self));
    memcpy(data, &self->size, sizeof(int));
    memcpy(data + sizeof(int), self->cost, self->size * self->size * sizeof(double));
    return data;
}

void Graph_fromBytes(Graph* self, void* bytes) {
    self->size = ((int*)bytes)[0];
    int doublesSize = self->size * self->size * sizeof(double);
    self->cost = (double*)malloc(doublesSize);
    memcpy(self->cost, bytes + sizeof(int), doublesSize);
}

void DFS(List* tree, int vertex, bool* visited, List* path) {
    List* self = tree + vertex;
    visited[vertex] = true;
    List_add(path, vertex);
    for (int i = 0; i < self->last; ++i) {
        int neighbor = self->data[i];
        if (!visited[neighbor]) {
            DFS(tree, neighbor, visited, path);
        }
    }
}

SalesmanPath Graph_getSimpleSalesmanPath(Graph* self) {
    Edge* edges = malloc(self->size * self->size * sizeof(Edge));
    DisjontSet* sets = malloc(self->size * sizeof(DisjontSet));

    bool* visited = malloc(self->size * sizeof(int));

    int index = 0;
    for (int i = 0; i < self-> size; ++i) {
        DisjontSet_init(&sets[i], i);
        visited[i] = false;
        for (int j = i + 1; j < self-> size; ++j) {
             edges[index++] = Edge_create(i, j, Graph_cost(self, i, j));
        }
    }

    qsort(edges, index, sizeof(Edge), Edge_compare);

    int merged = 1;
    index = 0;

    List* tree = malloc(self->size * sizeof(List));
    List path;
    List_init(&path, self->size);
    for (int i = 0; i < self->size; ++i) {
        List_init(tree + i, self->size);
    }

    while (merged < self->size) {        
        Edge edge = edges[index++];
        if (DisjontSet_representatn(sets + edge.a) != DisjontSet_representatn(sets + edge.b)) {
            DisjontSet_union(sets + edge.a, sets + edge.b);
            List_add(tree + edge.a, edge.b);
            List_add(tree + edge.b, edge.a);
            ++merged;
        }
    }

    DFS(tree, 0, visited, &path);
    SalesmanPath result;
    result.size = self->size;
    result.cost = 0;
    result.vertexes = malloc(self->size * sizeof(int));
    int prev = -1, current;
    for (int i = 0; i < self->size; ++i) {
        current = path.data[i];
        result.vertexes[i] = current;
        if (prev != -1) {
            result.cost += Graph_cost(self, prev, current);
        } else {
            result.cost += Graph_cost(self, current, path.data[self->size-1]);
        }
        prev = current;
    }

    List_free(&path);
    for (int i = 0; i < self->size; ++i) {
        List_free(tree + i);
    }
    free(tree);
    free(sets);
    free(edges);
    free(visited);

    return result;
}


