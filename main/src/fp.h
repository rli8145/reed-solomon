#pragma once

#include <gmpxx.h>
#include <stdexcept>

// Element of prime field F_p
class Fp {
public:
    mpz_class val;
    mpz_class p;

    Fp() : val(0), p(0) {}
    Fp(const mpz_class& value, const mpz_class& prime);

    Fp operator+(const Fp& rhs) const;
    Fp operator-(const Fp& rhs) const;
    Fp operator*(const Fp& rhs) const;
    Fp operator/(const Fp& rhs) const;
    Fp operator-() const;
    bool operator==(const Fp& rhs) const;
    bool operator!=(const Fp& rhs) const;

    // Modular inverse via extended Euclidean algorithm
    Fp inv() const;

    // Legendre symbol: returns -1, 0, or 1
    int legendre() const;

    // Square root via Tonelli-Shanks, throws if val is not a QR
    Fp sqrt() const;
};

// Modular exponentiation: base^exp (square-and-multiply)
Fp pow(const Fp& base, mpz_class exp);
