#pragma once

#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {
    void PromoteToHighOrder(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
    void ProjectMesh(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
}