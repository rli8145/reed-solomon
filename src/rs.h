#pragma once

#include "fp.h"
#include "poly.h"
#include <vector>
#include <optional>

class RSCode {
public:
    mpz_class p;              // field prime
    size_t k;                 // message length
    size_t n;                 // codeword length (n > k, n <= p)
    size_t t;                 // error correction capacity = floor((n-k)/2)
    std::vector<Fp> alphas;   // evaluation points

    // Default evaluation points: 0, 1, ..., n-1
    RSCode(const mpz_class& prime, size_t k, size_t n);

    // Custom evaluation points
    RSCode(const mpz_class& prime, size_t k, size_t n, std::vector<Fp> alphas);

    // Encode: interpret message as polynomial coefficients, evaluate at all n points
    std::vector<Fp> encode(const std::vector<Fp>& message) const;

    // Berlekamp-Welch decoding - returns nullopt if decoding fails
    std::optional<std::vector<Fp>> decode_bw(const std::vector<Fp>& received) const;
};
