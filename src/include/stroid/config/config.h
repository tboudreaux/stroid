#pragma once

#include <cstdint>

namespace stroid::config {

    struct OptimizationMethods {
        std::optional<bool> tmop{false};
        std::optional<bool> smoothstep{true};
    };

    /**
     * @brief Configuration parameters for stroid mesh generation.
     *
     * These values are typically loaded via
     * `fourdst::config::Config<stroid::config::MeshConfig>` from a TOML file.
     * The README shows the expected TOML layout under the `[main]` table.
     * Unspecified keys use the defaults defined here.
     */
    struct MeshConfig {
        /**
         * @brief Number of uniform refinement passes applied after topology creation.
         * @section toml
         * - [main].refinement_levels
         */
        std::optional<int> refinement_levels = 4;
        /**
         * @brief Polynomial order for high-order elements.
         * @section toml
         * - [main].order
         */
        std::optional<int> order = 3;
        /**
         * @brief Whether to include an external domain extending to `r_infinity`.
         * @section toml
         * - [main].include_external_domain
         */
        std::optional<bool> include_external_domain = true;

        /**
         * @brief Radius of the stellar core region.
         * @section toml
         * - [main].r_core
         */
        std::optional<double> r_core = 0.25;
        /**
         * @brief Radius of the stellar surface.
         * @section toml
         * - [main].r_star
         */
        std::optional<double> r_star = 1.0;
        /**
         * @brief Flattening factor for spheroidal shaping (0 = spherical, >0 = oblate).
         * @section toml
         * - [main].flattening
         */
        std::optional<double> flattening = 0;

        /**
         * @brief Outer radius of the external domain when enabled.
         * @section toml
         * - [main].r_infinity
         */
        std::optional<double> r_infinity = 6.0;

        /**
         * @brief Radius inside which transformations are skipped to avoid singularities.
         * @section toml
         * - [main].r_instability
         */
        std::optional<double> r_instability = 1e-14;
        /**
         * @brief Controls the smoothness/steepness of the core-to-envelope transition.
         * @section toml
         * - [main].core_steepness
         */
        std::optional<double> core_steepness = 1.0;

        /**
         * @brief Boundary attribute id for stellar surface
         * @section toml
         * - [main].surface_bdr_id
         */
        std::optional<size_t> surface_bdr_id = 1;

        /**
         * @brief Boundary attribute id for infinity in kelvin mapping
         * @section toml
         * - [main].inf_bdr_id
         */
        std::optional<size_t> inf_bdr_id = 2;

        /**
         * @brief Material attribute id for the core region
         * @section toml
         * - [main].core_id
         */
        std::optional<size_t> core_id = 1;

        /**
         * @brief Material attribute id for the envelope region
         * @section toml
         * - [main].envelope_id
         */
        std::optional<size_t> envelope_id = 2;

        /**
         * @brief Material attribute id for the external domain (if enabled)
         * @section toml
         * - [main].vacuum_id
         */
        std::optional<size_t> vacuum_id = 3;

        std::optional<OptimizationMethods> optimization_methods = OptimizationMethods{true, true};
    };
}
