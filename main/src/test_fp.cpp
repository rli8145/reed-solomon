#include "fp.h"
#include <gtest/gtest.h>

// Some primes for testing
static const mpz_class P17(17);
static const mpz_class P97(97);
static const mpz_class P7919(7919);
// 2^127 - 1 (Mersenne prime)
static const mpz_class P_BIG("170141183460469231731687303715884105727");

// Basic construction
TEST(FpConstruction, ReducesModP) {
    Fp a(20, P17);
    EXPECT_EQ(a.val, 3);  // 20 mod 17 = 3
}

TEST(FpConstruction, NegativeValueWraps) {
    Fp a(-1, P17);
    EXPECT_EQ(a.val, 16);  // -1 mod 17 = 16
}

TEST(FpConstruction, ZeroStaysZero) {
    Fp a(0, P17);
    EXPECT_EQ(a.val, 0);
}

// Arithmetic
TEST(FpArithmetic, Addition) {
    Fp a(10, P17), b(12, P17);
    EXPECT_EQ((a + b).val, 5);  // 22 mod 17 = 5
}

TEST(FpArithmetic, Subtraction) {
    Fp a(3, P17), b(10, P17);
    EXPECT_EQ((a - b).val, 10);  // 3 - 10 = -7 ≡ 10 (mod 17)
}

TEST(FpArithmetic, Multiplication) {
    Fp a(5, P17), b(7, P17);
    EXPECT_EQ((a * b).val, 1);  // 35 mod 17 = 1
}

TEST(FpArithmetic, UnaryNegation) {
    Fp a(5, P17);
    Fp neg = -a;
    EXPECT_EQ(neg.val, 12);  // -5 ≡ 12 (mod 17)
    EXPECT_EQ((a + neg).val, 0);
}

TEST(FpArithmetic, UnaryNegationOfZero) {
    Fp z(0, P17);
    EXPECT_EQ((-z).val, 0);
}

// Division and inverse
TEST(FpInverse, BasicInverse) {
    Fp a(3, P17);
    Fp inv = a.inv();
    EXPECT_EQ((a * inv).val, 1);
}

TEST(FpInverse, InverseOfZeroThrows) {
    Fp z(0, P17);
    EXPECT_THROW(z.inv(), std::domain_error);
}

TEST(FpInverse, AllNonzeroElementsInvertible) {
    for (int i = 1; i < 17; i++) {
        Fp a(i, P17);
        Fp inv = a.inv();
        EXPECT_EQ((a * inv).val, 1) << "Failed for i=" << i;
    }
}

TEST(FpDivision, DivisionIsMultByInverse) {
    Fp a(7, P97), b(13, P97);
    Fp quotient = a / b;
    EXPECT_EQ((quotient * b).val, a.val);
}


// Exponentiation
TEST(FpPow, ZeroExponent) {
    Fp a(42, P97);
    EXPECT_EQ(pow(a, 0).val, 1);
}

TEST(FpPow, FirstPower) {
    Fp a(42, P97);
    EXPECT_EQ(pow(a, 1).val, a.val);
}

TEST(FpPow, FermatsLittleTheorem) {
    // a^p ≡ a (mod p) for all a
    for (int i = 0; i < 97; i++) {
        Fp a(i, P97);
        EXPECT_EQ(pow(a, P97).val, a.val) << "Failed for i=" << i;
    }
}

TEST(FpPow, FermatInverse) {
    // a^(p-1) = 1 for a != 0
    Fp a(42, P7919);
    EXPECT_EQ(pow(a, P7919 - 1).val, 1);
}

TEST(FpPow, NegativeExponent) {
    Fp a(5, P17);
    Fp result = pow(a, -1);
    EXPECT_EQ((a * result).val, 1);
}

TEST(FpPow, BigPrime) {
    Fp a(123456789, P_BIG);
    EXPECT_EQ(pow(a, P_BIG - 1).val, 1);  // Fermat
}

// Legendre symbol
TEST(FpLegendre, ZeroIsZero) {
    Fp z(0, P17);
    EXPECT_EQ(z.legendre(), 0);
}

TEST(FpLegendre, OneIsQR) {
    Fp one(1, P17);
    EXPECT_EQ(one.legendre(), 1);
}

TEST(FpLegendre, QuadraticResiduesMod17) {
    // QRs mod 17: {1, 2, 4, 8, 9, 13, 15, 16}
    int qr[] = {1, 2, 4, 8, 9, 13, 15, 16};
    for (int v : qr) {
        Fp a(v, P17);
        EXPECT_EQ(a.legendre(), 1) << "Expected QR for " << v;
    }
}

TEST(FpLegendre, NonResiduesMod17) {
    int nr[] = {3, 5, 6, 7, 10, 11, 12, 14};
    for (int v : nr) {
        Fp a(v, P17);
        EXPECT_EQ(a.legendre(), -1) << "Expected NR for " << v;
    }
}

// Square root (Tonelli-Shanks)
TEST(FpSqrt, SqrtOfZero) {
    Fp z(0, P17);
    EXPECT_EQ(z.sqrt().val, 0);
}

TEST(FpSqrt, SqrtOfOne) {
    Fp one(1, P17);
    Fp r = one.sqrt();
    EXPECT_EQ((r * r).val, 1);
}

TEST(FpSqrt, AllQRsMod17) {
    int qr[] = {1, 2, 4, 8, 9, 13, 15, 16};
    for (int v : qr) {
        Fp a(v, P17);
        Fp r = a.sqrt();
        EXPECT_EQ((r * r).val, a.val) << "sqrt(" << v << ")^2 != " << v;
    }
}

TEST(FpSqrt, NonResidueThrows) {
    Fp a(3, P17);  // 3 is a NR mod 17
    EXPECT_THROW(a.sqrt(), std::domain_error);
}

TEST(FpSqrt, LargePrimeP3Mod4) {
    // P_BIG ≡ 3 (mod 4), so the fast path is taken.
    // Use a known square: 49 = 7^2.
    Fp a(49, P_BIG);
    EXPECT_EQ(a.legendre(), 1);
    Fp r = a.sqrt();
    EXPECT_EQ((r * r).val, a.val);
}

TEST(FpSqrt, PrimeP1Mod4) {
    // p = 7919 ≡ 3 (mod 4) actually, so pick p = 41 (≡ 1 mod 4)
    mpz_class p41(41);
    // QRs mod 41 include 2 (since 17^2 = 289 = 7*41 + 2)
    Fp a(2, p41);
    EXPECT_EQ(a.legendre(), 1);
    Fp r = a.sqrt();
    EXPECT_EQ((r * r).val, a.val);
}

TEST(FpSqrt, FullSweepMod41) {
    // 41 ≡ 1 (mod 4), exercises the full Tonelli-Shanks path.
    mpz_class p41(41);
    for (int v = 1; v < 41; v++) {
        Fp a(v, p41);
        if (a.legendre() == 1) {
            Fp r = a.sqrt();
            EXPECT_EQ((r * r).val, a.val)
                << "sqrt(" << v << ")^2 != " << v << " (mod 41)";
        }
    }
}

// Edge cases
TEST(FpEdge, PMinusOne) {
    Fp a(16, P17);  // p-1
    EXPECT_EQ((-a).val, 1);
    EXPECT_EQ((a * a).val, 1);  // (-1)^2 = 1
}

TEST(FpEdge, EqualityAcrossReduction) {
    Fp a(20, P17);
    Fp b(3, P17);
    EXPECT_EQ(a, b);  // both should be 3 mod 17
}

TEST(FpEdge, MismatchedPrimesThrow) {
    Fp a(1, P17);
    Fp b(1, P97);
    EXPECT_THROW(a + b, std::invalid_argument);
}
