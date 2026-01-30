#pragma once

namespace stroid::config {
    struct MeshConfig {
        int refinement_levels = 3;
        int order = 2;
        bool include_external_domain = false;

        double r_core = 2.5;
        double r_star = 5.0;
        double flattening = 0;

        double r_infinity = 6.0;

        double r_instability = 1e-14;
        double core_steepness = 1.0;
    };
}
