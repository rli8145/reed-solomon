#pragma once

#include "fp.h"
#include <vector>
#include <utility>
#include <tuple>

// Polynomial over Fp[x]
// Stored as coefficients where coeffs[i] is the coefficient of x^i (zero polynomial has an empty coefficient vector)

class Poly {
public:
    std::vector<Fp> coeffs;
    mpz_class p;  // the prime for the coefficient field

    Poly() : p(0) {}
    Poly(const mpz_class& prime) : p(prime) {}
    Poly(std::vector<Fp> coefficients, const mpz_class& prime);

    int degree() const;

    bool is_zero() const;

    // Leading coefficient (throws if zero polynomial)
    Fp leading() const;

    // Evaluate at a point using Horner's method
    Fp evaluate(const Fp& x) const;

    // Make polynomial monic (leading coefficient = 1)
    Poly monic() const;

    Poly operator+(const Poly& rhs) const;
    Poly operator-(const Poly& rhs) const;
    Poly operator*(const Poly& rhs) const;
    Poly operator*(const Fp& scalar) const;
    bool operator==(const Poly& rhs) const;
    bool operator!=(const Poly& rhs) const;

private:
    // Remove trailing zero coeffs
    void strip();
};

// Polynomial division: returns (quotient, remainder) s.t. a = b*q + r, deg r < deg b
std::pair<Poly, Poly> divmod(const Poly& a, const Poly& b);

// GCD via the Euclidean algorithm (result is monic)
Poly gcd(const Poly& a, const Poly& b);

// Extended GCD: returns (g, s, t) s.t. g = s*a + t*b (g is monic)
std::tuple<Poly, Poly, Poly> extended_gcd(const Poly& a, const Poly& b);

// Lagrange interpolation: given n points (x_i, y_i), returns unique
// polynomial of degree < n passing through all of them
Poly lagrange_interpolate(const std::vector<std::pair<Fp, Fp>>& points);
