#pragma once

#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {
    void ApplyEquiangular(mfem::Vector& pos);

    void ApplySpheroidal(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config);

    void ApplyKelvin(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config);

    void TransformPoint(mfem::Vector& pos, const fourdst::config::Config<config::MeshConfig> &config, int attribute_id);
}