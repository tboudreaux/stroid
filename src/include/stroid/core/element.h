#pragma once
#include <array>
#include <cstdint>

namespace stroid::core {
    struct HexElement {
        std::array<uint64_t, 8> vertices;
        int attribute_id;
    };
}