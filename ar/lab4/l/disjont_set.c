typedef struct __DisjontSet {
    struct __DisjontSet* parent;
    int id;
} DisjontSet;

void DisjontSet_init(DisjontSet* self, int id) {
    self->parent = NULL;
    self->id = id;
}

void* DisjontSet_representatn(DisjontSet* self) {
    if (self->parent == NULL) {
        return self;
    } else {
        self->parent = DisjontSet_representatn(self->parent);
        return self->parent;
    }
}

void DisjontSet_union(DisjontSet* a, DisjontSet* b) {
    DisjontSet* aRoot = DisjontSet_representatn(a);
    DisjontSet* bRoot = DisjontSet_representatn(b);
    aRoot->parent = bRoot;

}
