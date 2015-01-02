typedef struct __List {
    int last;
    int* data;
} List;

void List_init(List* self, int size) {
    self->last = 0;
    self->data = malloc(size * sizeof(int));
}

void List_add(List* self, int value) {
    self->data[self->last++] = value;
}

void List_free(List* self) {
    free(self->data);
}
