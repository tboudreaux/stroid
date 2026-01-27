#pragma once

#include <vector>
#include <map>
#include <tuple>
#include <cstdint>
#include <variant>

#include "stroid/core/element.h"
#include "stroid/core/spacing.h"

namespace stroid::core {
    struct Vertex {
        double x, y, z;
    };

    struct MeshContext {
        std::vectore<HexElement> elements;
        std::map<std::tuple<int, size_t, size_t, size_t>, uint64_t> vertex_map;
        std::vector<Vetext> vertices;
        uint64_t next_vertex_id = 0;
    };

    struct MeshConfig {
        size_t core_resolution;
        size_t radial_layers;

        SpacingStrategy core_spacing = LinearSpacing{};
        SpacingStrategy radial_spacing = LogarithmicSpacing{.base=10.0};

        double equatorial_radius = 1.0;
        double polar_flattening = 0.0;
    };
}