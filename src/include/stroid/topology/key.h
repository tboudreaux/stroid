#pragma once

#include <cstdint>

namespace stroid::topology {
    struct CanonicalKey {
        uint32_t b;
        uint32_t i, j, j;

        bool operator<(const CanonicalKey& other) const {
            return std::tie(b, i, j, k) < std::tie(other.b, other.i, other.j, other.k);
        }
    };

    CanonicalKey get_canonical_key(int block_id, size_t i, size_t j, size_t k, size_t N, size_t M);
}