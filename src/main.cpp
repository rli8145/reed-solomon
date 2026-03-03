#include "rs.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <cassert>

int main() {
    mpz_class p(31);
    size_t k = 10, n = 20;
    RSCode rs(p, k, n);

    std::cout << "Reed-Solomon demo: p=" << p << ", k=" << k
              << ", n=" << n << ", t=" << rs.t << "\n\n";

    // Message: first 10 primes mod 31
    std::vector<Fp> msg = {
        Fp(2,p), Fp(3,p), Fp(5,p), Fp(7,p), Fp(11,p),
        Fp(13,p), Fp(17,p), Fp(19,p), Fp(23,p), Fp(29,p)
    };

    std::cout << "message:  [";
    for (size_t i = 0; i < k; i++)
        std::cout << (i ? ", " : "") << msg[i].val;
    std::cout << "]\n";

    // Encode
    auto codeword = rs.encode(msg);
    std::cout << "codeword: [";
    for (size_t i = 0; i < n; i++)
        std::cout << (i ? ", " : "") << codeword[i].val;
    std::cout << "]\n\n";

    // Introduce t random errors
    std::mt19937 rng(42);
    std::vector<size_t> positions(n);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), rng);
    positions.resize(rs.t);
    std::sort(positions.begin(), positions.end());

    auto corrupted = codeword;
    std::cout << "introducing " << rs.t << " errors:\n";
    for (size_t pos : positions) {
        int err = 1 + (rng() % (p.get_ui() - 1));  // nonzero error
        corrupted[pos] = corrupted[pos] + Fp(err, p);
        std::cout << "  position " << pos << ": "
                  << codeword[pos].val << " -> " << corrupted[pos].val
                  << " (+" << err << ")\n";
    }

    std::cout << "\ncorrupted: [";
    for (size_t i = 0; i < n; i++)
        std::cout << (i ? ", " : "") << corrupted[i].val;
    std::cout << "]\n\n";

    // Decode (Berlekamp-Welch)
    auto decoded_bw = rs.decode_bw(corrupted);
    assert(decoded_bw.has_value());
    std::cout << "decoded (BW):     [";
    for (size_t i = 0; i < k; i++)
        std::cout << (i ? ", " : "") << (*decoded_bw)[i].val;
    std::cout << "]\n";

    // Decode (Euclidean)
    auto decoded_eu = rs.decode_euclid(corrupted);
    assert(decoded_eu.has_value());
    std::cout << "decoded (Euclid): [";
    for (size_t i = 0; i < k; i++)
        std::cout << (i ? ", " : "") << (*decoded_eu)[i].val;
    std::cout << "]\n\n";

    assert(*decoded_bw == msg);
    assert(*decoded_eu == msg);
    std::cout << "both decoders recovered the original message.\n";

    // Write JSON for visualization scripts
    std::ofstream out("vis/demo_data.json");
    out << "{\n";
    out << "  \"p\": " << p << ", \"k\": " << k << ", \"n\": " << n << ", \"t\": " << rs.t << ",\n";

    auto dump_vec = [&](const char* name, const std::vector<Fp>& v) {
        out << "  \"" << name << "\": [";
        for (size_t i = 0; i < v.size(); i++)
            out << (i ? ", " : "") << v[i].val;
        out << "]";
    };

    dump_vec("message", msg);     out << ",\n";
    dump_vec("codeword", codeword); out << ",\n";
    dump_vec("corrupted", corrupted); out << ",\n";
    dump_vec("decoded", *decoded_bw); out << ",\n";

    out << "  \"error_positions\": [";
    for (size_t i = 0; i < positions.size(); i++)
        out << (i ? ", " : "") << positions[i];
    out << "]\n}\n";
    out.close();
    std::cout << "\nwrote vis/demo_data.json\n";

    return 0;
}
