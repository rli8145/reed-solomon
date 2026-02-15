#include "fp.h"

// helpers 

// Reduce v into [0, p).
static mpz_class mod(const mpz_class& v, const mpz_class& p) {
    mpz_class r = v % p;
    if (r < 0) r += p;
    return r;
}

static void check_field(const Fp& a, const Fp& b) {
    if (a.p != b.p)
        throw std::invalid_argument("Fp: mismatched primes");
}

// construction
Fp::Fp(const mpz_class& value, const mpz_class& prime)
    : val(mod(value, prime)), p(prime) {}

// arithmetic

Fp Fp::operator+(const Fp& rhs) const {
    check_field(*this, rhs);
    return Fp(val + rhs.val, p);
}

Fp Fp::operator-(const Fp& rhs) const {
    check_field(*this, rhs);
    return Fp(val - rhs.val, p);
}

Fp Fp::operator*(const Fp& rhs) const {
    check_field(*this, rhs);
    return Fp(val * rhs.val, p);
}

Fp Fp::operator/(const Fp& rhs) const {
    check_field(*this, rhs);
    return *this * rhs.inv();
}

Fp Fp::operator-() const {
    return Fp(p - val, p);
}

bool Fp::operator==(const Fp& rhs) const {
    return p == rhs.p && val == rhs.val;
}

bool Fp::operator!=(const Fp& rhs) const {
    return !(*this == rhs);
}

// modular inverse (extended GCD)

Fp Fp::inv() const {
    if (val == 0)
        throw std::domain_error("Fp: inverse of zero");

    mpz_class g, x, y;
    mpz_gcdext(g.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(),
               val.get_mpz_t(), p.get_mpz_t());
    if (g != 1)
        throw std::domain_error("Fp: element not invertible");
    return Fp(x, p);
}

// exponentiation (square-and-multiply)

Fp pow(const Fp& base, mpz_class exp) {
    if (base.p == 0)
        throw std::invalid_argument("Fp: uninitialized element");

    // handle negative exponents
    if (exp < 0) {
        return pow(base.inv(), -exp);
    }

    Fp result(1, base.p);
    Fp b = base;

    while (exp > 0) {
        if (mpz_odd_p(exp.get_mpz_t()))
            result = result * b;
        b = b * b;
        exp >>= 1;
    }
    return result;
}

// Legendre symbol

int Fp::legendre() const {
    // (a/p) = a^((p-1)/2) mod p   →   0, 1, or p-1 (≡ -1)
    Fp r = pow(*this, (p - 1) / 2);
    if (r.val == 0) return 0;
    if (r.val == 1) return 1;
    return -1;
}

// Tonelli-Shanks square root

Fp Fp::sqrt() const {
    if (val == 0) return *this;

    int leg = legendre();
    if (leg == -1)
        throw std::domain_error("Fp::sqrt: not a quadratic residue");

    // Special case: p ≡ 3 (mod 4)  →  sqrt = a^((p+1)/4)
    if ((p % 4) == 3)
        return pow(*this, (p + 1) / 4);

    // Factor p-1 = Q * 2^S with Q odd
    mpz_class Q = p - 1;
    unsigned long S = 0;
    while (mpz_even_p(Q.get_mpz_t())) {
        Q >>= 1;
        S++;
    }

    // Find a quadratic non-residue z
    Fp z(2, p);
    while (z.legendre() != -1)
        z = Fp(z.val + 1, p);

    Fp M_unused(S, p);          // only need S as integer
    unsigned long M = S;
    Fp c = pow(z, Q);           // c = z^Q
    Fp t = pow(*this, Q);       // t = n^Q
    Fp R = pow(*this, (Q + 1) / 2);

    while (true) {
        if (t.val == 1) return R;

        // Find least i such that t^(2^i) = 1
        unsigned long i = 0;
        Fp tmp = t;
        while (tmp.val != 1) {
            tmp = tmp * tmp;
            i++;
        }

        // Update
        Fp b = c;
        for (unsigned long j = 0; j < M - i - 1; j++)
            b = b * b;

        M = i;
        c = b * b;
        t = t * c;
        R = R * b;
    }
}
