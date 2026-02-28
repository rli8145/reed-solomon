#include "fp.h"
#include "poly.h"
#include "gaussian_elim.h"
#include <cassert>
#include <iostream>

static void test_fp() {
    mpz_class p(17);

    // Arithmetic spot checks
    assert(Fp(10, p) + Fp(12, p) == Fp(5, p));   // 22 mod 17
    assert(Fp(3, p) - Fp(10, p) == Fp(10, p));   // -7 mod 17

    // Inverse correctness
    for (int i = 1; i < 17; i++) {
        Fp a(i, p);
        assert(a * a.inv() == Fp(1, p));
    }

    std::cout << "  fp: ok\n";
}

static void test_poly() {
    mpz_class p(7);

    // Division: a = b*q + r, deg r < deg b
    auto a = Poly({Fp(3,p), Fp(1,p), Fp(4,p), Fp(2,p)}, p);
    auto b = Poly({Fp(1,p), Fp(0,p), Fp(3,p)}, p);
    auto [q, r] = divmod(a, b);
    assert(b * q + r == a);
    assert(r.degree() < b.degree());

    // GCD of (x^3 - 1) and (x^2 - 1) over F7 is (x - 1)
    auto f = Poly({Fp(6,p), Fp(0,p), Fp(0,p), Fp(1,p)}, p);
    auto g = Poly({Fp(6,p), Fp(0,p), Fp(1,p)}, p);
    auto d = gcd(f, g);
    assert(d == Poly({Fp(6,p), Fp(1,p)}, p));  // x + 6 = x - 1

    // Lagrange round-trip: interpolate a known polynomial from its evaluations
    auto orig = Poly({Fp(1,p), Fp(6,p), Fp(3,p), Fp(2,p)}, p);
    std::vector<std::pair<Fp,Fp>> pts;
    for (int i = 0; i < 4; i++) {
        Fp x(i, p);
        pts.push_back({x, orig.evaluate(x)});
    }
    assert(lagrange_interpolate(pts) == orig);

    std::cout << "  poly: ok\n";
}

static void test_gaussian_elim() {
    mpz_class p(7);

    // System 1: 2x + 3y = 1, 4x + y = 5  (mod 7)
    // By hand: row2 -= 2*row1 => y-row gives (1 - 2*3)y = 5 - 2*1 => -5y = 3 => 2y = 3 => y = 3*4 = 12 = 5
    // Then 2x = 1 - 3*5 = 1 - 15 = -14 = 0 => x = 0
    // Check: 2*0 + 3*5 = 15 = 1 mod 7, 4*0 + 5 = 5. Correct.
    {
        std::vector<std::vector<Fp>> A = {
            {Fp(2,p), Fp(3,p)},
            {Fp(4,p), Fp(1,p)}
        };
        std::vector<Fp> b = {Fp(1,p), Fp(5,p)};
        auto x = solve(A, b);
        assert(x.has_value());
        assert((*x)[0] == Fp(0,p));
        assert((*x)[1] == Fp(5,p));
    }

    // System 2: 3x3 system
    // x + 2y + 3z = 1,  2x + 5y + 3z = 2,  x + 0y + 6z = 3  (mod 7)
    // Verify by substitution
    {
        std::vector<std::vector<Fp>> A = {
            {Fp(1,p), Fp(2,p), Fp(3,p)},
            {Fp(2,p), Fp(5,p), Fp(3,p)},
            {Fp(1,p), Fp(0,p), Fp(6,p)}
        };
        std::vector<Fp> b = {Fp(1,p), Fp(2,p), Fp(3,p)};
        auto x = solve(A, b);
        assert(x.has_value());
        // Check Ax = b
        for (size_t i = 0; i < 3; i++) {
            Fp sum(0, p);
            for (size_t j = 0; j < 3; j++)
                sum = sum + A[i][j] * (*x)[j];
            assert(sum == b[i]);
        }
    }

    // Singular system: two identical rows
    {
        std::vector<std::vector<Fp>> A = {
            {Fp(1,p), Fp(2,p)},
            {Fp(1,p), Fp(2,p)}
        };
        std::vector<Fp> b = {Fp(3,p), Fp(3,p)};
        auto x = solve(A, b);
        assert(!x.has_value());
    }

    std::cout << "  gaussian_elim: ok\n";
}

int main() {
    std::cout << "tests:\n";
    test_fp();
    test_poly();
    test_gaussian_elim();
    std::cout << "all passed.\n";
    return 0;
}
