#include "gaussian_elim.h"
#include <stdexcept>
#include <utility>

std::optional<std::vector<Fp>> solve(
    std::vector<std::vector<Fp>> A, std::vector<Fp> b) {
    size_t n = A.size();
    if (n == 0) return std::vector<Fp>{};
    if (b.size() != n)
        throw std::invalid_argument("solve: dimension mismatch");
    for (auto& row : A)
        if (row.size() != n)
            throw std::invalid_argument("solve: non-square matrix");

    const mpz_class& p = b[0].p;

    // Forward elimination
    for (size_t col = 0; col < n; col++) {
        // Partial pivot: find a row with nonzero entry in this column
        size_t pivot = n;
        for (size_t row = col; row < n; row++) {
            if (A[row][col].val != 0) { pivot = row; break; }
        }
        if (pivot == n) return std::nullopt;  // singular

        if (pivot != col) {
            std::swap(A[col], A[pivot]);
            std::swap(b[col], b[pivot]);
        }

        Fp inv = A[col][col].inv();

        // Eliminate below
        for (size_t row = col + 1; row < n; row++) {
            Fp factor = A[row][col] * inv;
            for (size_t j = col; j < n; j++)
                A[row][j] = A[row][j] - factor * A[col][j];
            b[row] = b[row] - factor * b[col];
        }
    }

    // Back substitution
    Vec x(n, Fp(0, p));
    for (int i = (int)n - 1; i >= 0; i--) {
        Fp sum = b[i];
        for (size_t j = i + 1; j < n; j++)
            sum = sum - A[i][j] * x[j];
        x[i] = sum * A[i][i].inv();
    }

    return x;
}
