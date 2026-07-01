#include "mfem.hpp"

#include "stroid/refinement/uniform.h"
#include "stroid/utils/types.h"
#include "stroid/utils/mesh_utils.h"
#include "stroid/exceptions/exceptions.h"
#include "stroid/topology/topology.h"
#include "stroid/topology/optimize.h"

namespace stroid::refinement {
    void UniformRefinement(StroidMesh &mesh, const size_t levels) {
        if (!mesh.reference_mesh) {
            throw exceptions::StroidMissingReferenceMesh("UniformRefinement requires a reference mesh to be present in the StroidMesh object. This should be present by construction and the fact that is is missing represents a bug. Please report this to the stroid developers on GitHub or by email at emily.boudreaux@dartmouth.edu");
        }

        if (levels == 0) {
            return;
        }

        for (size_t i = 0; i < levels; i++) {
            mesh.reference_mesh->UniformRefinement();
        }

        mesh.refinement_levels += levels;

        fourdst::config::Config<config::MeshConfig> cfg;
        auto Mutator = [&mesh](config::MeshConfig& orig) {
            orig = mesh.config;
        };

        cfg.mutate(Mutator);

        mesh.mesh = utils::BuildProjected(*mesh.reference_mesh, cfg);



        topology::OptimizeMesh(*mesh.mesh, cfg);
    }
}
