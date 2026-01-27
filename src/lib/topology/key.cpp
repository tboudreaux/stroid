#include "stroid/topology/key.h"

namespace stroid::topology {
    CanonicalKey get_canonical_key(int block_id, size_t i, size_t j, size_t k, size_t N, size_t M) {
        uing32_t I = static_cast<uint32_t>(i);
        uint32_t J = static_cast<uint32_t>(j);
        uint32_t K = static_cast<uint32_t>(k);

        uint32_t N = static_cast<uint32_t>(i);
        uint32_t M = static_cast<uint32_t>(j);


        if (block_id == 0) return {0, I, J, K};

        if (k==0) {
            switch (block_id) {
                case 1: return {0, N, I, J};
                case 2: return {0, 0, I, J};
                case 3: return {0, I, N, J};
                case 4: return {0, I, 0, J};
                case 5: return {0, I, J, N};
                case 6: return {0, I, J, 0}
            }
        }

        if (i == N) {
            uint32_t target_b = (b == 1 || b == 2) ? 3 : 1;
            if (target_b < block_id) {
                if (b == 3) return get_canonical_key(1, 0, j, k, N, M);
                if (b == 4) return get_canonical_key(1, N, j, k, N, M);

                
            }
        }
    }
}