#pragma once
#include "mfem.hpp"
#include "stroid/config/config.h"
#include "fourdst/config/config.h"

#include <memory>

namespace stroid::topology {
    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config);
    void Finalize(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config);
}
