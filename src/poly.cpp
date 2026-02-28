#include "poly.h"
#include <stdexcept>
#include <algorithm>

// helpers 

static void check_field(const Poly& a, const Poly& b) {
    if (a.p != b.p)
        throw std::invalid_argument("Poly: mismatched primes");
}

// construction 

Poly::Poly(std::vector<Fp> coefficients, const mpz_class& prime)
    : coeffs(std::move(coefficients)), p(prime) {
    strip();
}

void Poly::strip() {
    while (!coeffs.empty() && coeffs.back().val == 0)
        coeffs.pop_back();
}

// queries

int Poly::degree() const {
    return static_cast<int>(coeffs.size()) - 1;
}

bool Poly::is_zero() const {
    return coeffs.empty();
}

Fp Poly::leading() const {
    if (is_zero())
        throw std::domain_error("Poly: zero polynomial has no leading coefficient");
    return coeffs.back();
}

// evaluation (Horner's method)

Fp Poly::evaluate(const Fp& x) const {
    if (is_zero())
        return Fp(0, p);

    Fp result = coeffs.back();
    for (int i = degree() - 1; i >= 0; --i)
        result = result * x + coeffs[i];
    return result;
}

//  monic

Poly Poly::monic() const {
    if (is_zero())
        throw std::domain_error("Poly: zero polynomial cannot be made monic");
    Fp lc_inv = leading().inv();
    std::vector<Fp> c(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i)
        c[i] = coeffs[i] * lc_inv;
    return Poly(std::move(c), p);
}

// arithmetic

Poly Poly::operator+(const Poly& rhs) const {
    check_field(*this, rhs);
    size_t n = std::max(coeffs.size(), rhs.coeffs.size());
    std::vector<Fp> c(n, Fp(0, p));
    for (size_t i = 0; i < coeffs.size(); ++i)
        c[i] = c[i] + coeffs[i];
    for (size_t i = 0; i < rhs.coeffs.size(); ++i)
        c[i] = c[i] + rhs.coeffs[i];
    return Poly(std::move(c), p);
}

Poly Poly::operator-(const Poly& rhs) const {
    check_field(*this, rhs);
    size_t n = std::max(coeffs.size(), rhs.coeffs.size());
    std::vector<Fp> c(n, Fp(0, p));
    for (size_t i = 0; i < coeffs.size(); ++i)
        c[i] = c[i] + coeffs[i];
    for (size_t i = 0; i < rhs.coeffs.size(); ++i)
        c[i] = c[i] - rhs.coeffs[i];
    return Poly(std::move(c), p);
}

Poly Poly::operator*(const Poly& rhs) const {
    check_field(*this, rhs);
    if (is_zero() || rhs.is_zero())
        return Poly(p);

    size_t n = coeffs.size() + rhs.coeffs.size() - 1;
    std::vector<Fp> c(n, Fp(0, p));
    for (size_t i = 0; i < coeffs.size(); ++i)
        for (size_t j = 0; j < rhs.coeffs.size(); ++j)
            c[i + j] = c[i + j] + coeffs[i] * rhs.coeffs[j];
    return Poly(std::move(c), p);
}

Poly Poly::operator*(const Fp& scalar) const {
    if (scalar.val == 0)
        return Poly(p);
    std::vector<Fp> c(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i)
        c[i] = coeffs[i] * scalar;
    return Poly(std::move(c), p);
}

bool Poly::operator==(const Poly& rhs) const {
    if (p != rhs.p) return false;
    return coeffs.size() == rhs.coeffs.size() &&
           std::equal(coeffs.begin(), coeffs.end(), rhs.coeffs.begin());
}

bool Poly::operator!=(const Poly& rhs) const {
    return !(*this == rhs);
}

// divmod (polynomial long division)

std::pair<Poly, Poly> divmod(const Poly& a, const Poly& b) {
    if (b.is_zero())
        throw std::domain_error("Poly divmod: division by zero polynomial");
    check_field(a, b);

    if (a.degree() < b.degree())
        return {Poly(a.p), a};

    const mpz_class& p = a.p;
    Fp lc_inv = b.leading().inv();

    // Work with coefficient arrays
    std::vector<Fp> rem = a.coeffs;
    int deg_diff = a.degree() - b.degree();
    std::vector<Fp> quot(deg_diff + 1, Fp(0, p));

    for (int i = deg_diff; i >= 0; --i) {
        Fp coeff = rem[i + b.degree()] * lc_inv;
        quot[i] = coeff;
        for (int j = 0; j <= b.degree(); ++j)
            rem[i + j] = rem[i + j] - coeff * b.coeffs[j];
    }

    return {Poly(std::move(quot), p), Poly(std::move(rem), p)};
}

// GCD (Euclidean algorithm)
Poly gcd(const Poly& a, const Poly& b) {
    if (b.is_zero())
        return a.is_zero() ? a : a.monic();
    auto [q, r] = divmod(a, b);
    return gcd(b, r);
}

// Extended GCD

std::tuple<Poly, Poly, Poly> extended_gcd(const Poly& a, const Poly& b) {
    const mpz_class& p = a.p;

    if (b.is_zero()) {
        if (a.is_zero())
            return {Poly(p), Poly(p), Poly(p)};
        // g = a (made monic), s = 1/lc(a), t = 0
        Fp lc_inv = a.leading().inv();
        Poly s({lc_inv}, p);
        return {a.monic(), s, Poly(p)};
    }

    auto [q, r] = divmod(a, b);
    auto [g, s1, t1] = extended_gcd(b, r);
    // g = s1*b + t1*r = s1*b + t1*(a - q*b) = t1*a + (s1 - t1*q)*b
    Poly s = t1;
    Poly t = s1 - t1 * q;
    return {g, s, t};
}

// Lagrange interpolation

Poly lagrange_interpolate(const std::vector<std::pair<Fp, Fp>>& points) {
    if (points.empty())
        throw std::invalid_argument("Poly lagrange: empty point set");

    const mpz_class& p = points[0].first.p;
    size_t n = points.size();

    // Start with the zero polynomial
    Poly result(p);

    for (size_t i = 0; i < n; ++i) {
        // Build the Lagrange basis polynomial L_i(x)
        // L_i(x) = product_{j != i} (x - x_j) / (x_i - x_j)

        // Numerator: product of (x - x_j) for j != i
        Poly basis({Fp(1, p)}, p);  // start with constant 1
        Fp denom(1, p);

        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            // Multiply basis by (x - x_j)
            Poly linear({-points[j].first, Fp(1, p)}, p);
            basis = basis * linear;
            denom = denom * (points[i].first - points[j].first);
        }

        // L_i(x) = basis / denom, scaled by y_i
        Fp scale = points[i].second * denom.inv();
        result = result + basis * scale;
    }

    return result;
}
