#include <fstream>
#include <array>
#include "create_system.h"
#include "functions.h"
#include "solution_SLAE/method_Gauss_Seidel.h"

int main() {
    const std::size_t N = 150;
    const SLAE<double> slae = system<double>(N - 2, f, phi);
    const Vector<double> x_0{(N - 2) * (N - 2)};
    const double eps = 1e-11;

    const Vector<double> res = method_Gauss_Seidel(slae.A_, slae.b_, x_0, 10000, eps);

    const double step = 1. / (N - 1);
    std::ofstream data("data_for_heat_map.csv");
    data << N << '\n';

    for (std::size_t i = 0; i < N; ++i) {
        data << phi(i * step, 0) << " ";
    }
    data << '\n';
    for (std::size_t i = 1; i < N - 1; ++i) {
        data << phi(0, i * step) << " ";
        for (std::size_t j = 1; j < N - 1; ++j) {
            data << res((i - 1) * (N - 2) + (j - 1)) << " ";
        }
        data << phi((N - 1) * step, i * step) << "\n";
    }
    for (std::size_t i = 0; i < N; ++i) {
        data << phi(i * step, (N - 1) * step) << " ";
    }
}