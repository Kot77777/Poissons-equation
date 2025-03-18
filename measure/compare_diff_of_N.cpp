#include <fstream>
#include <array>
#include "create_system.h"
#include "functions.h"
#include "solution_SLAE/method_Gauss_Seidel.h"

template<typename T>
T get_err(const Vector<T> &vec, const std::size_t N){
    T large{};
    T diff{};
    const T step = 1. / (N + 2);
    for (std::size_t i = 1; i <= N; ++i) {
        for (std::size_t j = 1; j <= N; ++j) {
            diff = fabs(vec((i - 1) * N + (j - 1)) - sin(j * step));
            if (diff > large) {
                large = diff;
            }
        }
    }
    return large;

}

int main() {
    const std::array<std::size_t, 4> N{10, 20, 40, 80};
    const double eps = 1e-11;
    std::ofstream data("data.csv");
    data << "N" << "," << "err" << '\n';

    for (std::size_t i : N) {
        SLAE<double> slae = system<double>(i - 2, f, phi);
        Vector<double> x_0{(i - 2) * (i - 2)};
        Vector<double> res = method_Gauss_Seidel(slae.A_, slae.b_, x_0, 10000, eps);
        data << i << "," << get_err(res, i - 2) << "\n";
    }
}