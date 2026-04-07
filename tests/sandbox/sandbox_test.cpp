#include "fourdst/config/config.h"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"
#include "stroid/topology/topology.h"

#include "mfem.hpp"

#include <print>

#include "stroid/topology/optimize.h"

struct SandboxConfig {
    std::string host = "localhost";
    int port = 19916;
    bool visualize = true;
};

using MeshConfig = fourdst::config::Config<stroid::config::MeshConfig>;
using UserConfig = fourdst::config::Config<SandboxConfig>;

int main() {
    MeshConfig mesh_cfg;
    mesh_cfg.load("default.toml");

    const UserConfig user_cfg;


    std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(mesh_cfg);
    stroid::topology::Finalize(*mesh, mesh_cfg);
    stroid::topology::PromoteToHighOrder(*mesh, mesh_cfg);
    stroid::topology::ProjectMesh(*mesh, mesh_cfg);

    if (mesh_cfg->optimization_methods.has_value() && mesh_cfg->optimization_methods.value().tmop.has_value() && mesh_cfg->optimization_methods.value().tmop.value()) {
        stroid::topology::ApplyTMOP(*mesh, mesh_cfg);
    }

    stroid::IO::ViewMesh(*mesh, "Sandbox Mesh", stroid::IO::VISUALIZATION_MODE::ELEMENT_ID, user_cfg->host, user_cfg->port);
    stroid::IO::SaveMesh(*mesh, "sandbox.mesh");

     return 0;


}