#pragma once

#include "fp.h"
#include <vector>
#include <optional>

// Solve Ax = b over Fp via Gaussian elimination with partial pivoting.
// Returns nullopt if the system is singular.
std::optional<std::vector<Fp>> solve(
    std::vector<std::vector<Fp>> A, std::vector<Fp> b);
