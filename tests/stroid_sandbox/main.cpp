#include <print>
#include "stroid/topology/topology.h"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include <memory>
#include "mfem.hpp"
#include "stroid/topology/curvilinear.h"
#include "stroid/utils/mesh_utils.h"

#include "fourdst/config/config.h"

int main() {
    const fourdst::config::Config<stroid::config::MeshConfig> cfg;

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);
    //
    // stroid::utils::MarkFlippedElements(*mesh);
    // stroid::utils::MarkFlippedBoundaryElements(*mesh);


    stroid::IO::ViewMesh(*mesh, "Spheroidal Mesh", stroid::IO::VISUALIZATION_MODE::BOUNDARY_ELEMENT_ID);
    // stroid::IO::VisualizeFaceValence(*mesh);
    // stroid::IO::SaveVTU(*mesh, "SpheroidalMesh");
}
