#ifndef SLAE_H
#define SLAE_H
#include <primitives/CSR_matrix.h>
#include "primitives/vector_from_vector.h"

template<typename T>
struct SLAE {
    CSR_matrix<T> A_;
    Vector<T> b_;
};

#endif //SLAE_H
