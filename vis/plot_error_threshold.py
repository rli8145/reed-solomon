#!/usr/bin/env python3
"""Plot theoretical vs empirical error correction threshold.

For fixed p and k, varies n and plots:
  - Theoretical max correctable errors t = floor((n-k)/2) vs n
  - Empirical boundary: for each n, introduces increasing numbers of errors
    and records whether decoding succeeds (via the C++ demo binary).

Requires the C++ 'test_threshold' binary (built by this script).
"""

import subprocess
import json
import tempfile
import os
import matplotlib.pyplot as plt
import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Write a small C++ program that tests decoding for given parameters
THRESHOLD_SRC = r"""
#include "rs.h"
#include <iostream>
#include <random>
#include <algorithm>

int main(int argc, char** argv) {
    // Usage: test_threshold <p> <k> <n_min> <n_max> <trials>
    mpz_class p(argv[1]);
    size_t k = std::stoul(argv[2]);
    size_t n_min = std::stoul(argv[3]);
    size_t n_max = std::stoul(argv[4]);
    int trials = std::stoi(argv[5]);

    std::mt19937 rng(123);

    // JSON output
    std::cout << "[\n";
    bool first_outer = true;

    for (size_t n = n_min; n <= n_max; n++) {
        RSCode rs(p, k, n);
        size_t t = rs.t;

        // Generate a random message
        std::vector<Fp> msg;
        for (size_t i = 0; i < k; i++)
            msg.emplace_back(rng() % p.get_ui(), p);
        auto codeword = rs.encode(msg);

        if (!first_outer) std::cout << ",\n";
        first_outer = false;
        std::cout << "  {\"n\": " << n << ", \"t\": " << t << ", \"results\": [";

        bool first_inner = true;
        // Test 0 through n-k errors
        for (size_t num_err = 0; num_err <= n - k; num_err++) {
            int successes = 0;
            for (int trial = 0; trial < trials; trial++) {
                auto corrupted = codeword;

                // Pick num_err random positions
                std::vector<size_t> pos(n);
                std::iota(pos.begin(), pos.end(), 0);
                std::shuffle(pos.begin(), pos.end(), rng);
                pos.resize(num_err);

                for (size_t pi : pos) {
                    int err = 1 + (rng() % (p.get_ui() - 1));
                    corrupted[pi] = corrupted[pi] + Fp(err, p);
                }

                auto decoded = rs.decode_euclid(corrupted);
                if (decoded.has_value() && *decoded == msg)
                    successes++;
            }
            if (!first_inner) std::cout << ", ";
            first_inner = false;
            std::cout << successes;
        }
        std::cout << "]}";
    }
    std::cout << "\n]\n";
    return 0;
}
"""

def build_threshold_binary():
    """Write and compile the threshold test binary."""
    src_path = os.path.join(PROJ, "src", "test_threshold.cpp")
    with open(src_path, "w") as f:
        f.write(THRESHOLD_SRC)

    # Build using cmake
    subprocess.run(
        ["cmake", "--build", os.path.join(PROJ, "build"), "--target", "test_threshold"],
        cwd=PROJ, check=True, capture_output=True
    )
    os.remove(src_path)
    return os.path.join(PROJ, "build", "test_threshold")


def main():
    p = 31
    k = 10
    n_min = k + 2   # need n > k, and at least 2 redundancy for t >= 1
    n_max = 28       # must be <= p
    trials = 20

    # Add temporary build target
    cmake_path = os.path.join(PROJ, "CMakeLists.txt")
    with open(cmake_path) as f:
        cmake_orig = f.read()

    cmake_extra = "\nadd_executable(test_threshold src/test_threshold.cpp)\ntarget_link_libraries(test_threshold PRIVATE rs_lib)\n"

    # Write threshold source and update CMakeLists
    src_path = os.path.join(PROJ, "src", "test_threshold.cpp")
    with open(src_path, "w") as f:
        f.write(THRESHOLD_SRC)
    with open(cmake_path, "a") as f:
        f.write(cmake_extra)

    try:
        subprocess.run(
            ["cmake", "-S", PROJ, "-B", os.path.join(PROJ, "build"),
             "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"],
            cwd=PROJ, check=True, capture_output=True
        )
        subprocess.run(
            ["cmake", "--build", os.path.join(PROJ, "build"), "--target", "test_threshold"],
            cwd=PROJ, check=True, capture_output=True
        )

        binary = os.path.join(PROJ, "build", "test_threshold")
        result = subprocess.run(
            [binary, str(p), str(k), str(n_min), str(n_max), str(trials)],
            capture_output=True, text=True, check=True
        )
        data = json.loads(result.stdout)
    finally:
        # Restore CMakeLists and clean up
        with open(cmake_path, "w") as f:
            f.write(cmake_orig)
        if os.path.exists(src_path):
            os.remove(src_path)
        # Reconfigure to remove stale target
        subprocess.run(
            ["cmake", "-S", PROJ, "-B", os.path.join(PROJ, "build"),
             "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"],
            cwd=PROJ, capture_output=True
        )

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    ns = [d["n"] for d in data]
    ts = [d["t"] for d in data]

    # Theoretical threshold
    ax.plot(ns, ts, "k-o", linewidth=2, markersize=6, label="theoretical t = ⌊(n−k)/2⌋", zorder=5)

    # Empirical heatmap: for each n, plot success rate vs number of errors
    for d in data:
        n_val = d["n"]
        results = d["results"]
        for num_err, successes in enumerate(results):
            rate = successes / trials
            color = plt.cm.RdYlGn(rate)
            ax.plot(n_val, num_err, "s", color=color, markersize=8, alpha=0.8)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("decoding success rate ({} trials)".format(trials))

    ax.set_xlabel("codeword length n")
    ax.set_ylabel("number of errors")
    ax.set_title("RS error correction: p={}, k={}\nblack line = theoretical limit".format(p, k))
    ax.legend(loc="upper left")

    plt.tight_layout()
    plt.savefig("vis/error_threshold.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
