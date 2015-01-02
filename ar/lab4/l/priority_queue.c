typedef struct __QNode {
    double priority;
    void * value;
} QNode;

typedef struct __PriorityQueue {
    QNode* nodes;
    int size;
    int availableSize;
} PriorityQueue;

void PriorityQueue_init(PriorityQueue& self) {
    QNode* = malloc()
}

void PriorityQueue_free(PriorityQueue& self) {
    free(self->nodes);
}

void PriorityQueue_insert(PriorityQueue& self, QNode& value) {

}

QNode& PriorityQueue_pop(PriorityQueue& self) {

}


