#pragma once

#include "mfem.hpp"

#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::utils {
    /**
     * @brief Mark elements with negative Jacobian determinant.
     * @details Elements detected as flipped are assigned attribute 999.
     * @param mesh Mesh to scan and update in-place.
     */
    void MarkFlippedElements(mfem::Mesh& mesh);
    /**
     * @brief Mark boundary elements whose outward normal orientation is flipped.
     * @details Boundary elements detected as flipped are assigned attribute 500.
     * @param mesh Mesh to scan and update in-place.
     */
    void MarkFlippedBoundaryElements(mfem::Mesh& mesh);

    void ExportJacobianRadialProfile(mfem::Mesh& mesh, const std::string& filename);

    std::unique_ptr<mfem::Mesh> BuildProjected(const mfem::Mesh& reference, const fourdst::config::Config<config::MeshConfig>& cfg);

}