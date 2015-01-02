#include "string.h"
#include "stdbool.h"

int sign(double value) {
    return (value > 0) - (value < 0);
}

long perm_compute_all(int len, int preflen) {
    long result = 1, mn = len;
    for(int i = 0; i < preflen; ++i) {
        result *= mn;
        --mn;
    }
    return result;
}

void set_perm(int * elems, int len, long ord, int preflen) {
    int* copy = malloc(len * sizeof(int));
    bool* available = malloc(len * sizeof(bool));
    long all, index, count;
    int lenCopy = len;

    memcpy(copy, elems, len * sizeof(int));
    memset(available, true, len);

    for (int i = 0; i < lenCopy; ++i) {
        --preflen, --len;
        if (preflen >= 0) {
            all = perm_compute_all(len, preflen);
            index = ord / all;
            ord %= all;
        } else {
            index = 0;
        }
        count = 0;
        for (int j = 0; j < lenCopy; ++j) {
            if (available[j] && ++count > index) {
                available[j] = false;
                elems[i] = copy[j];
                break;
            }
        }
    }

    free(copy);
    free(available);
}

void identity(int* tab, int size) {
    for(int i = 0; i < size; ++i) {
        tab[i] = i;
    }
}

void random_perm(int* tab, int size) {
    for(int i = 1, j; i < size; ++i) {
        j = rand() % i;
        tab[i]^=tab[j];
        tab[j]^=tab[i];
        tab[i]^=tab[j];
    }
}
