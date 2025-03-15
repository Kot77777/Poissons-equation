#include <gtest/gtest.h>
#include "create_system.h"
#include "primitives/CSR_matrix.h"
#include "solution_SLAE/method_Gauss_Seidel.h"

const auto phi = [](const double x, const double y) {
  return sin(x);
};

const auto f = [](const double x, const double y) {
    return -sin(x);
};

TEST(create_system, system) {
    const std::size_t N = 5;
    const SLAE<double> slae = system<double>(N - 2, f, phi);
    const CSR_matrix<double> A{slae.A_};
    const Vector<double> b{slae.b_};
    const Vector<double> x_0{(N - 2) * (N - 2)};
    const double eps = 10e-15;

    const Vector<double> res = method_Gauss_Seidel(A, b, x_0, 1000, eps);

    for (std::size_t i = 0; i < (N - 2) * (N - 2); ++i) {
        for (std::size_t j = 0; j < (N - 2) * (N - 2); ++j) {
            std::cout << A(i, j) << " ";
        }
        std::cout << '\n';
    }

    for(std::size_t i = 0; i < (N - 2) * (N - 2); ++i) {
        std::cout << res(i) << '\n';
    }
}
