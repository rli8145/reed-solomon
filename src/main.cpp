#include "rs.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>

int main() {
    // p = 131: smallest prime >= 128, so every ASCII byte fits in F_p
    mpz_class p(131);

    std::cout << "reed-solomon demo (F_" << p << ")\n";
    std::cout << "enter a message (max 65 chars): ";

    std::string input;
    if (!std::getline(std::cin, input) || input.empty()) {
        std::cerr << "no input\n";
        return 1;
    }

    size_t k = input.size();
    if (k > 65) {
        std::cerr << "message too long (max 65 chars, need room for redundancy with n <= p=131)\n";
        return 1;
    }

    // n = min(2k, 130) gives ~50% redundancy, capped so n <= p
    size_t n = std::min(2 * k, (size_t)130);
    if (n <= k) n = k + 2;
    RSCode rs(p, k, n);

    std::cout << "\n--- parameters ---\n";
    std::cout << "  p = " << p << "  (field size)\n";
    std::cout << "  k = " << k << "  (message symbols)\n";
    std::cout << "  n = " << n << "  (codeword length)\n";
    std::cout << "  t = " << rs.t << "  (correctable errors)\n";

    // Convert string to Fp vector
    std::vector<Fp> msg;
    msg.reserve(k);
    for (char c : input)
        msg.emplace_back((unsigned char)c, p);

    std::cout << "\n--- message ---\n";
    std::cout << "  string:  \"" << input << "\"\n";
    std::cout << "  values:  [";
    for (size_t i = 0; i < k; i++)
        std::cout << (i ? ", " : "") << msg[i].val;
    std::cout << "]\n";

    // Encode
    auto codeword = rs.encode(msg);
    std::cout << "\n--- codeword (" << n << " symbols) ---\n  [";
    for (size_t i = 0; i < n; i++)
        std::cout << (i ? ", " : "") << codeword[i].val;
    std::cout << "]\n";

    // Introduce t random errors
    std::random_device rd;
    std::mt19937 rng(rd());

    auto corrupted = codeword;
    std::vector<size_t> positions(n);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), rng);
    positions.resize(rs.t);
    std::sort(positions.begin(), positions.end());

    std::cout << "\n--- introducing " << rs.t << " errors ---\n";
    for (size_t pos : positions) {
        int err = 1 + (rng() % (p.get_ui() - 1));
        corrupted[pos] = corrupted[pos] + Fp(err, p);
        std::cout << "  position " << pos << ": "
                  << codeword[pos].val << " -> " << corrupted[pos].val << "\n";
    }

    // Decode
    std::cout << "\n--- decoding ---\n";

    auto decoded_bw = rs.decode_bw(corrupted);
    auto decoded_eu = rs.decode_euclid(corrupted);

    auto to_string = [&](const std::vector<Fp>& v) {
        std::string s;
        for (auto& f : v) s += (char)f.val.get_ui();
        return s;
    };

    if (decoded_bw.has_value()) {
        std::string recovered = to_string(*decoded_bw);
        std::cout << "  berlekamp-welch: \"" << recovered << "\""
                  << (*decoded_bw == msg ? "  [correct]" : "  [WRONG]") << "\n";
    } else {
        std::cout << "  berlekamp-welch: failed\n";
    }

    if (decoded_eu.has_value()) {
        std::string recovered = to_string(*decoded_eu);
        std::cout << "  euclidean:       \"" << recovered << "\""
                  << (*decoded_eu == msg ? "  [correct]" : "  [WRONG]") << "\n";
    } else {
        std::cout << "  euclidean:       failed\n";
    }

    // Write JSON for visualization script
    std::ofstream out("vis/demo_data.json");
    auto dump_vec = [&](const char* name, const std::vector<Fp>& v) {
        out << "  \"" << name << "\": [";
        for (size_t i = 0; i < v.size(); i++)
            out << (i ? ", " : "") << v[i].val;
        out << "]";
    };
    out << "{\n";
    out << "  \"message_str\": \"" << input << "\",\n";
    out << "  \"p\": " << p << ", \"k\": " << k << ", \"n\": " << n << ", \"t\": " << rs.t << ",\n";
    dump_vec("codeword", codeword); out << ",\n";
    dump_vec("corrupted", corrupted); out << ",\n";
    out << "  \"error_positions\": [";
    for (size_t i = 0; i < positions.size(); i++)
        out << (i ? ", " : "") << positions[i];
    out << "]\n}\n";
    std::cout << "\nwrote vis/demo_data.json\n";

    return 0;
}
