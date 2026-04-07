#pragma once

#include "mfem.hpp"
#include "fourdst/config/base.h"
#include "stroid/config/config.h"


namespace stroid::topology {

    /**
     * @breif Apply target matrix optimization to improve conditioning of the mesh
     */
    void ApplyTMOP(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
}
