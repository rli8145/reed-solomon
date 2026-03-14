# Reed-Solomon Error Correcting Codes

A from-scratch C++ implementation of Reed-Solomon error correcting codes over finite fields 𝔽ₚ, built on top of GMP for exact arithmetic. Includes two independent decoders - Berlekamp-Welch and Euclidean - and Python visualizations via matplotlib.

## Components

A longer writeup can be found [here](reed-solomon.pdf) (shoutout to MATH145 and MATH146).

### `fp.h / fp.cpp` — Finite field 𝔽ₚ

Elements of ℤ/pℤ for prime p, backed by GMP's `mpz_class`. All arithmetic reduces mod p. Multiplicative inverse via the extended Euclidean algorithm. Also provides Legendre symbol, Tonelli-Shanks square root, and modular exponentiation.

### `poly.h / poly.cpp` — Polynomials over 𝔽ₚ[x]

Polynomials with coefficients in 𝔽ₚ, stored as a coefficient vector (`coeffs[i]` = coefficient of xⁱ).

| Function | Description |
|---|---|
| `evaluate(x)` | Horner's method |
| `divmod(a, b)` | Polynomial long division; returns (q, r) with a = b·q + r, deg r < deg b |
| `gcd(a, b)` | Euclidean algorithm; result is monic |
| `extended_gcd(a, b)` | Returns (g, s, t) with g = s·a + t·b (Bézout's identity) |
| `lagrange_interpolate(pts)` | Unique polynomial of degree < n through n given points |

### `gaussian_elim.h / gaussian_elim.cpp` — Linear algebra over 𝔽ₚ

`solve(A, b)` solves an n×n linear system Ax = b over 𝔽ₚ using Gaussian elimination with partial pivoting. Division is exact (no floating point), so partial pivoting is only needed to avoid dividing by zero. Returns `nullopt` if the system is singular.

### `rs.h / rs.cpp` — Reed-Solomon codes

`RSCode(p, k, n)` constructs a code over 𝔽ₚ with message length k, codeword length n, and evaluation points α₀, …, αₙ₋₁ = 0, …, n-1. The error correction capacity is t = ⌊(n-k)/2⌋. Custom evaluation points can be supplied via the second constructor.

**Encoding.** The message (m₀, …, mₖ₋₁) is interpreted as coefficients of a polynomial m(x) of degree < k. The codeword is (m(α₀), …, m(αₙ₋₁)).

**Berlekamp-Welch decoder** (`decode_bw`). Finds polynomials E(x) (error locator, monic, degree e) and Q(x) = m(x)·E(x) (degree < k+e) satisfying yᵢ·E(αᵢ) = Q(αᵢ) at all n received positions. This gives a linear system in the unknown coefficients of E and Q. The decoder tries e = t, t-1, …, 0 until the system is non-singular, then recovers m = Q/E.

**Euclidean decoder** (`decode_euclid`). Lagrange-interpolates the received values to get a polynomial R(x) of degree < n, then builds the vanishing polynomial V(x) = ∏ᵢ(x − αᵢ). Runs the extended Euclidean algorithm on (V, R), stopping as soon as deg(remainder) < (n+k)/2. At that point the remainder is proportional to Q = m·E and the corresponding Bézout coefficient is proportional to E, so m = remainder / coefficient.

Both decoders return `nullopt` if decoding fails i.e. more than t errors.

## Dependencies

- C++17
- [GMP](https://gmplib.org/) (`brew install gmp`)
- CMake ≥ 3.16
- Python3 with numpy, matplotlib (`pip install -r requirements.txt`)

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Usage

### Interactive demo

Encodes a user-supplied ASCII string as a Reed-Solomon codeword, introduces the maximum number of correctable errors, and recovers the original message with both decoders.

```bash
./build/demo
# enter a message (max 65 chars): hello world
```

The demo uses p = 131 (smallest prime ≥ 128) so each ASCII character maps directly to a field element. For a message of length k, the codeword has length n = min(2k, 130), giving error correction capacity t = ⌊(n-k)/2⌋.

Sample output for the message `"Hello! My name is Ryan"`:

```
reed-solomon demo (F_131)
enter a message (max 65 chars): Hello! My name is Ryan

--- parameters ---
  p = 131  (field size)
  k = 22  (message symbols)
  n = 44  (codeword length)
  t = 11  (correctable errors)

--- message ---
  string:  "Hello! My name is Ryan"
  values:  [72, 101, 108, 108, 111, 33, 32, 77, 121, 32, 110, 97, 109, 101, 32, 105, 115, 32, 82, 121, 97, 110]

--- codeword (44 symbols) ---
  [72, 72, 26, 73, 75, 117, 45, 114, 22, 8, 122, 94, 97, 18, 106, 90, 76, 129, 90, 105, 26, 37, 70, 50, 43, 97, 31, 28, 84, 98, 124, 17, 118, 110, 4, 87, 9, 97, 42, 128, 44, 75, 37, 91]

--- introducing 11 errors ---
  position 0: 72 -> 59
  ...

--- decoding ---
  berlekamp-welch: "Hello! My name is Ryan"  [correct]
  euclidean:       "Hello! My name is Ryan"  [correct]
```

After running the demo, generate the bar chart visualization:

```bash
python3 vis/demo_recovery.py    # saves vis/demo_recovery.png
```

### Error threshold plot

Empirically measures how many errors can be corrected across varying codeword lengths, compared against the theoretical limit.

```bash
python3 vis/plot_error_threshold.py    # saves vis/error_threshold.png
```

### Tests

```bash
./build/test
```

Covers: 𝔽ₚ arithmetic and inverses, polynomial division/GCD/Lagrange interpolation, Gaussian elimination, and RS encode/decode with 0 errors, t errors, and t+1 errors (expected failure).

---
