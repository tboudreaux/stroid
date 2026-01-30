#include "mfem.hpp"
#include <vector>
#include <memory>

#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {

    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config) {
        int nVert = config->include_external_domain ? 24 : 16;
        int nElem = config->include_external_domain ? 13 : 7;
        int nBev  = 6;

        auto mesh = std::make_unique<mfem::Mesh>(3, nVert, nElem, nBev, 3);

        auto add_box = [&](double scale) {
            for (const double z : {-scale, scale})
                for (const double y : {-scale, scale})
                    for (const double x : {-scale, scale})
                        mesh->AddVertex(x, y, z);
        };

        add_box(config->r_core);
        add_box(config->r_star);
        if (config->include_external_domain) {
            add_box(config->r_infinity);
        }

        const int core_v[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        mesh->AddHex(core_v, 1);

        int shells[6][8] = {
            {8, 9, 11, 10, 0, 1, 3, 2},
            {4, 5, 7, 6, 12, 13, 15, 14}, // +Z face
            {0, 1, 5, 4, 8, 9, 13, 12},   // -Y face
            {10, 11, 15, 14, 2, 3, 7, 6},
            {1, 3, 7, 5, 9, 11, 15, 13},  // +X face
            {0, 4, 6, 2, 8, 12, 14, 10}   // -X face
        };
        for (const auto & shell : shells) mesh->AddHex(shell, 2);

        const int bdr_quads[6][4] = {
            {12, 13, 15, 14},
            {13, 9, 11, 15},
            {9, 8, 10, 11},
            {8, 12, 14, 10},
            {8, 9, 13, 12},
            {14, 15, 11, 10}
        };

        for (const auto& bdr: bdr_quads) {
            mesh->AddBdrQuad(bdr, 1);
        }

        return mesh;
    }

    void Finalize(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        mesh.FinalizeTopology();
        mesh.Finalize();
        mesh.CheckElementOrientation(true);
        mesh.CheckBdrElementOrientation(true);
        for (int i = 0; i < config->refinement_levels; ++i) {
            mesh.UniformRefinement();
        }

        if (!mesh.Conforming()) {
            std::cerr << "WARNING: Mesh has been detected to be non conforming!" << std::endl;
        }


    }

}
