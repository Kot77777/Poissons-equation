#ifndef CREATE_SYSTEM_H
#define CREATE_SYSTEM_H
#include "SLAE.h"

template<typename T, typename Callable_1, typename Callable_2>
SLAE<T> system(const std::size_t N, const Callable_1 &f0, const Callable_2 &phi0) {

    const T h = 1. / (N + 2);
    const T h_kv = h * h;

    const auto phi = [&phi0, h](const std::size_t j, const std::size_t i) {
        return phi0(i * h, j * h);
    };

    const auto f = [&f0, h](const std::size_t j, const std::size_t i) {
      return f0(i * h, j * h);
    };

    Vector<T> b{N * N};
    std::map<std::array<std::size_t, 2>, T> DOK;


    for (std::size_t i = 0; i < N * N; ++i) {
        DOK[{i, i}] = -4. / h_kv;
    }

    for (std::size_t i = 0; i < N * N - (N - 1); i += N) {
        for (std::size_t j = 0; j < N - 1; ++j) {
            DOK[{i + j, i + j + 1}] = 1. / h_kv;
            DOK[{i + j + 1, i + j}] = 1. / h_kv;
        }
    }

    for (std::size_t i = 0; i < N * N - N; ++i) {
        DOK[{i, i + N}] = 1. / h_kv;
        DOK[{i + N, i}] = 1. / h_kv;
    }

    b(0) = f(1, 1) - (phi(0, 1) + phi(1, 0)) / h_kv;
    b(N - 1) = f(1, N) - (phi(0, N) + phi(1, N + 1)) / h_kv;
    b(N * N - N) = f(N, 1) - (phi(N, 0) + phi(N + 1, 1)) / h_kv;
    b(N * N - 1) = f(N, N) - (phi(N + 1, N) + phi(N, N + 1)) / h_kv;

    for (std::size_t i = 1; i < N - 1; ++i) {
        b(i) = f(1, i + 1) - phi(0, i + 1) / h_kv;
        b(i + N * (N - 1)) = f(N, i + 1) - phi(N + 1, i + 1) / h_kv;
        b(i * N) = f(i + 1, 1) - phi(i + 1, 0) / h_kv;
        b(N * i + (N - 1)) = f(i + 1, N) - phi(i + 1, N + 1) / h_kv;
    }

    for (std::size_t i = 1; i < N - 1; ++i) {
        for (std::size_t j = 1; j < N - 1; ++j) {
            b(N * i + j) = f(i + 1, j + 1);
        }
    }
    return SLAE<T>{CSR_matrix<T>{DOK, N*N, N*N}, b};
}

#endif //CREATE_SYSTEM_H
