#include "stroid/topology/mapping.h"
#include <cmath>
#include <algorithm>

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
        pos(2) *= (1.0 - config->flattening);
    }

    void ApplyKelvin(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config) {
        const double r = pos.Norml2();
        if (r <= config->r_star) {
            return;
        }

        double xi = (r - config->r_star) / (config->r_infinity - config->r_star);
        xi = std::min(0.999, std::max(0.0, xi)); // Clamp xi to [0, 0.999]

        const double r_new = config->r_star + xi / (1.0 - xi);
        const double scale = r_new / r;
        pos *= scale;
    }

    void TransformPoint(mfem::Vector &pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id) {
        double l_inf = 0.0;
        for (int i = 0; i < pos.Size(); ++i) {
            l_inf = std::max(l_inf, std::abs(pos(i)));
        }

        if (l_inf < config->r_instability) return;

        // Gnomonic projection
        const double r_log = pos.Norml2();
        mfem::Vector unit_dir = pos;
        unit_dir /= r_log;

        ApplyEquiangular(unit_dir);
        unit_dir /= unit_dir.Norml2(); // Re-normalize

        if (l_inf <= config->r_core) {
            const double t = l_inf / config->r_core;
            double alpha = std::pow(t, config->core_steepness);

            // Smoothstep function to apply C1 continuity
            alpha = alpha * alpha * (3.0 - 2.0 * alpha);

            mfem::Vector pos_cartesian = pos;
            mfem::Vector pos_spherical = unit_dir;

            pos_spherical *= l_inf;


            for (int d = 0; d < pos.Size(); ++d) {
                pos(d) = (1.0 - alpha) * pos_cartesian(d) + alpha * pos_spherical(d);
            }

            ApplySpheroidal(pos, config);
            return;
        }



        if (l_inf <= config->r_star) {
            const double xi = (l_inf - config->r_core) / (config->r_star - config->r_core);
            const double r_phys = config->r_core + xi * (config->r_star - config->r_core);

            pos = unit_dir;
            pos *= r_phys;

            ApplySpheroidal(pos, config);
        } else {
            pos = unit_dir;
            pos *= l_inf;

            ApplyKelvin(pos, config);
            ApplySpheroidal(pos, config);
        }
    }
}
