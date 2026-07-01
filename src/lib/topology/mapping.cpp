#include "stroid/topology/mapping.h"
#include "stroid/exceptions/exceptions.h"
#include <cmath>
#include <algorithm>
#include <array>
#include <utility>
#include <format>
#include <string>

namespace {
    template<int n, int k>
    consteval int nCr() {
        if constexpr (k > n) {
            return 0;
        } else {
            if constexpr (constexpr int kk = (k * 2 > n) ? (n - k) : k; kk == 0) {
                return 1;
            } else {
                int result = n;

                for (int i = 2; i <= kk; ++i) {
                    result *= (n - i + 1);
                    result /= i;
                }

                return result;
            }
        }
    }

    template <int n>
    double GeneralizedSmoothstep(const double x) {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;

        double sum = 0.0;

        auto compute_term = [&]<std::size_t k>(std::integral_constant<std::size_t, k>) {
            return nCr<n + k, k>() * std::pow(1.0 - x, k);
        };

        auto unroller = [&]<std::size_t... ks>(std::index_sequence<ks...>) {
            return (compute_term(std::integral_constant<std::size_t, ks>{}) + ...);
        };

        sum = unroller(std::make_index_sequence<n + 1>{});
        return sum * std::pow(x, n + 1);
    }

    template <std::size_t... Is>
    constexpr auto make_smoothstep_dispatch_table(std::index_sequence<Is...>) {
        return std::array<double(*)(double), sizeof...(Is)>{
            &GeneralizedSmoothstep<Is + 1>...
        };
    }

    constexpr int MAX_SMOOTHSTEP_ORDER = 10;
    constexpr auto smoothstep_dispatch = make_smoothstep_dispatch_table(
        std::make_index_sequence<MAX_SMOOTHSTEP_ORDER>{}
    );
}

namespace stroid::topology {
    void ApplyEquiangular(mfem::Vector &pos) {
        const double x = pos(0);
        const double y = pos(1);
        const double z = pos(2);

        const double absX = std::abs(x);
        const double absY = std::abs(y);
        const double absZ = std::abs(z);

        const double maxAbs = std::max({absX, absY, absZ});

        if (maxAbs < 1e-14) return;

        if (absX == maxAbs) {
            pos(1) = x * std::tan(M_PI / 4.0 * (y/x));
            pos(2) = x * std::tan(M_PI / 4.0 * (z/x));
        } else if (absY == maxAbs) {
            pos(0) = y * std::tan(M_PI / 4.0 * (x/y));
            pos(2) = y * std::tan(M_PI / 4.0 * (z/y));
        } else { // absZ == maxAbs
            pos(0) = z * std::tan(M_PI / 4.0 * (x/z));
            pos(1) = z * std::tan(M_PI / 4.0 * (y/z));
        }
    }

    void ApplySpheroidal(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config) {
        pos(2) *= (1.0 - config->flattening.value());
    }

    void TransformPoint(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id) {
        double X = pos(0);
        double Y = pos(1);
        double Z = pos(2);

        double maxAbs = std::max({std::abs(X), std::abs(Y), std::abs(Z)});
        if (maxAbs < 1e-14) return;

        double cx = X / maxAbs;
        double cy = Y / maxAbs;
        double cz = Z / maxAbs;

        double sx = cx * std::sqrt(1.0 - cy*cy/2.0 - cz*cz/2.0 + cy*cy*cz*cz/3.0);
        double sy = cy * std::sqrt(1.0 - cx*cx/2.0 - cz*cz/2.0 + cx*cx*cz*cz/3.0);
        double sz = cz * std::sqrt(1.0 - cx*cx/2.0 - cy*cy/2.0 + cx*cx*cy*cy/3.0);

        mfem::Vector unit_dir(3);
        unit_dir(0) = sx;
        unit_dir(1) = sy;
        unit_dir(2) = sz;

        if (maxAbs <= config->r_core.value()) {
            double nx = X / config->r_core.value();
            double ny = Y / config->r_core.value();
            double nz = Z / config->r_core.value();

            pos(0) = nx * std::sqrt(1.0 - ny*ny/2.0 - nz*nz/2.0 + ny*ny*nz*nz/3.0);
            pos(1) = ny * std::sqrt(1.0 - nx*nx/2.0 - nz*nz/2.0 + nx*nx*nz*nz/3.0);
            pos(2) = nz * std::sqrt(1.0 - nx*nx/2.0 - ny*ny/2.0 + nx*nx*ny*ny/3.0);

            pos *= config->r_core.value();

            ApplySpheroidal(pos, config);
            return;
        }

        if (maxAbs <= config->r_star.value()) {
            const double xi = (maxAbs - config->r_core.value()) / (config->r_star.value() - config->r_core.value());
            const double r_phys = config->r_core.value() + xi * (config->r_star.value() - config->r_core.value());

            pos = unit_dir;
            pos *= r_phys;

            ApplySpheroidal(pos, config);
        } else {
            pos = unit_dir;
            pos *= maxAbs;

            ApplySpheroidal(pos, config);
        }
    }

    // void TransformPoint(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id) {
    //     double l_inf = 0.0;
    //     for (int i = 0; i < pos.Size(); ++i) {
    //         l_inf = std::max(l_inf, std::abs(pos(i)));
    //     }
    //
    //     if (l_inf < config->r_instability) return;
    //
    //     // Gnomonic projection
    //     const double r_log = pos.Norml2();
    //     mfem::Vector unit_dir = pos;
    //     unit_dir /= r_log;
    //
    //     ApplyEquiangular(unit_dir);
    //     unit_dir /= unit_dir.Norml2(); // Re-normalize
    //
    //     if (l_inf <= config->r_core) {
    //         const double t = l_inf / config->r_core.value();
    //         double alpha = std::pow(t, config->core_steepness.value());
    //         const size_t order = config->continuity_order.value_or(2);
    //         if (order < 1 || order > MAX_SMOOTHSTEP_ORDER) {
    //             const std::string err_msg = std::format("Invalid continuity order: {}. Continuity order must be between (inclusive) 1 and {}. To push to higher orders you must update MAX_SMOOTHSTEP_ORDER in src/lib/topology/mapping.cpp and recompile.", order, MAX_SMOOTHSTEP_ORDER);
    //             throw exceptions::StroidContinuityError(err_msg);
    //         }
    //
    //         alpha = smoothstep_dispatch[order - 1](alpha); // We use this funky method as it keeps smoothstep calculation largely offloaded to compile time rather than run-time
    //
    //         mfem::Vector pos_cartesian = pos;
    //         mfem::Vector pos_spherical = unit_dir;
    //
    //         pos_spherical *= l_inf;
    //         bool run_smoothstep = false;
    //
    //
    //         if (config->optimization_methods.has_value() && config->optimization_methods.value().smoothstep.has_value() && config->optimization_methods.value().smoothstep.value()) {
    //             run_smoothstep = true;
    //         }
    //
    //
    //         if (run_smoothstep) {
    //             for (int d = 0; d < pos.Size(); ++d) {
    //                 pos(d) = (1.0 - alpha) * pos_cartesian(d) + alpha * pos_spherical(d);
    //             }
    //         }
    //
    //         ApplySpheroidal(pos, config);
    //         return;
    //     }
    //
    //     if (l_inf <= config->r_star) {
    //         const double xi = (l_inf - config->r_core.value()) / (config->r_star.value() - config->r_core.value());
    //         const double r_phys = config->r_core.value() + xi * (config->r_star.value() - config->r_core.value());
    //
    //         pos = unit_dir;
    //         pos *= r_phys;
    //
    //         ApplySpheroidal(pos, config);
    //     } else {
    //         pos = unit_dir;
    //         pos *= l_inf;
    //
    //         ApplySpheroidal(pos, config);
    //     }
    // }
}
