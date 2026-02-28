#include "rs.h"
#include "gaussian_elim.h"
#include <stdexcept>

RSCode::RSCode(const mpz_class& prime, size_t k, size_t n)
    : p(prime), k(k), n(n), t((n - k) / 2) {
    if (n <= k)
        throw std::invalid_argument("RSCode: need n > k");
    if (n > prime)
        throw std::invalid_argument("RSCode: need n <= p");
    alphas.reserve(n);
    for (size_t i = 0; i < n; i++)
        alphas.emplace_back(i, p);
}

RSCode::RSCode(const mpz_class& prime, size_t k, size_t n, std::vector<Fp> alphas)
    : p(prime), k(k), n(n), t((n - k) / 2), alphas(std::move(alphas)) {
    if (n <= k)
        throw std::invalid_argument("RSCode: need n > k");
    if (this->alphas.size() != n)
        throw std::invalid_argument("RSCode: need exactly n evaluation points");
}

std::vector<Fp> RSCode::encode(const std::vector<Fp>& message) const {
    if (message.size() != k)
        throw std::invalid_argument("RSCode::encode: message length must be k");

    Poly m(std::vector<Fp>(message.begin(), message.end()), p);
    std::vector<Fp> codeword(n);
    for (size_t i = 0; i < n; i++)
        codeword[i] = m.evaluate(alphas[i]);
    return codeword;
}

std::optional<std::vector<Fp>> RSCode::decode_bw(const std::vector<Fp>& received) const {
    if (received.size() != n)
        throw std::invalid_argument("RSCode::decode_bw: received length must be n");

    // Berlekamp-Welch: find E(x) monic of degree e and Q(x) of degree < k+e
    // such that for each i: yᵢ · E(αᵢ) = Q(αᵢ)
    //
    // We don't know how many errors occurred, so try e = t, t-1, ..., 0.
    // The system has k+2e unknowns (e coeffs of E + k+e coeffs of Q).
    // We use k+2e of the n equations. When e matches the actual error count,
    // the system is non-singular (BW guarantee).
    //
    // For each point i:
    //   yᵢ·(e₀ + e₁αᵢ + ... + e_{e-1}αᵢ^{e-1} + αᵢ^e) = q₀ + q₁αᵢ + ... + q_{k+e-1}αᵢ^{k+e-1}
    // Rearranged (moving the monic term to the RHS):
    //   e₀·yᵢ + ... + e_{e-1}·yᵢ·αᵢ^{e-1} - q₀ - q₁·αᵢ - ... = -yᵢ·αᵢ^e

    for (int e = (int)t; e >= 0; e--) {
        size_t num_unknowns = k + 2 * (size_t)e;

        std::vector<std::vector<Fp>> A(num_unknowns, std::vector<Fp>(num_unknowns, Fp(0, p)));
        std::vector<Fp> b(num_unknowns, Fp(0, p));

        for (size_t i = 0; i < num_unknowns; i++) {
            Fp yi = received[i];
            Fp alpha_pow(1, p);

            for (int j = 0; j < e; j++) {
                A[i][j] = yi * alpha_pow;
                alpha_pow = alpha_pow * alphas[i];
            }
            // alpha_pow = α_i^e
            b[i] = Fp(0, p) - yi * alpha_pow;

            Fp alpha_pow2(1, p);
            for (size_t j = 0; j < k + (size_t)e; j++) {
                A[i][e + j] = Fp(0, p) - alpha_pow2;
                alpha_pow2 = alpha_pow2 * alphas[i];
            }
        }

        auto sol = solve(A, b);
        if (!sol.has_value()) continue;

        // Extract E(x): e_0,...,e_{e-1}, plus leading 1
        std::vector<Fp> e_coeffs(e + 1);
        for (int j = 0; j < e; j++)
            e_coeffs[j] = (*sol)[j];
        e_coeffs[e] = Fp(1, p);
        Poly E(std::move(e_coeffs), p);

        // Extract Q(x)
        std::vector<Fp> q_coeffs(k + (size_t)e);
        for (size_t j = 0; j < k + (size_t)e; j++)
            q_coeffs[j] = (*sol)[e + j];
        Poly Q(std::move(q_coeffs), p);

        // Recover m(x) = Q(x) / E(x) (should divide exactly)
        if (E.is_zero()) continue;
        auto [m, rem] = divmod(Q, E);
        if (!rem.is_zero()) continue;
        if (m.degree() >= (int)k) continue;

        std::vector<Fp> msg(k, Fp(0, p));
        for (size_t i = 0; i < m.coeffs.size() && i < k; i++)
            msg[i] = m.coeffs[i];
        return msg;
    }

    return std::nullopt;
}
