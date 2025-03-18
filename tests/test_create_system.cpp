#include <gtest/gtest.h>
#include "create_system.h"
#include "primitives/CSR_matrix.h"
#include "solution_SLAE/method_Gauss_Seidel.h"
#include "functions.h"

TEST(create_system, system) {
    const std::size_t N = 4;
    const SLAE<double> slae = system<double>(N - 2, f, phi);
    const CSR_matrix<double> A{slae.A_};
    const Vector<double> b{slae.b_};
    const Vector<double> x_0{(N - 2) * (N - 2)};
    const double eps = 1e-12;

    const Vector<double> res = method_Gauss_Seidel(A, b, x_0, 10000, eps);

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
