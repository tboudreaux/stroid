#pragma once

#include <cstdint>
#include <vector>
#include "stroid/core/element.h"

namespace stroid::topology {
    class BlockTopology {
        BlockTopology(size_t nx, size_t ny, size_t nz);

        uint64_t get_vertex_id(size_t i, size_t j, size_t k) const;
        std::vector<core::HexElement> generate_elements(int attribute_id) const;
    };
}