#include "mfem.hpp"
#include <vector>
#include <memory>

#include "stroid/config/config.h"
#include "fourdst/config/config.h"

namespace stroid::topology {

    std::unique_ptr<mfem::Mesh> BuildSkeleton(const fourdst::config::Config<config::MeshConfig> & config) {
        int nVert = config->include_external_domain ? 24 : 16;
        int nElem = config->include_external_domain ? 13 : 7;
        int nBev  = config->include_external_domain ? 12 : 6;

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
        mesh->AddHex(core_v, config->core_id);

        std::vector<std::array<int, 8>> stellar_shells = {
            {8, 9, 11, 10, 0, 1, 3, 2},
            {4, 5, 7, 6, 12, 13, 15, 14}, // +Z face
            {0, 1, 5, 4, 8, 9, 13, 12},   // -Y face
            {10, 11, 15, 14, 2, 3, 7, 6},
            {1, 3, 7, 5, 9, 11, 15, 13},  // +X face
            {0, 4, 6, 2, 8, 12, 14, 10}   // -X face
        };
        for (const auto & shell : stellar_shells) {
            mesh->AddHex(shell.data(), config->envelope_id);
        }

        if (config->include_external_domain) {
            std::vector<std::array<int, 8>> vacuum_shells;
            vacuum_shells.push_back({8, 9, 13, 12, 16, 17, 21, 20});
            vacuum_shells.push_back({9, 11, 15, 13, 17, 19, 23, 21});
            vacuum_shells.push_back({11, 10, 14, 15, 19, 18, 22, 23});
            vacuum_shells.push_back({10, 8, 12, 14, 18, 16, 20, 22});
            vacuum_shells.push_back({12, 13, 15, 14, 20, 21, 23, 22});
            vacuum_shells.push_back({10, 11, 9, 8, 18, 19, 17, 16});
            for (const auto & shell : vacuum_shells) {
                mesh->AddHex(shell.data(), config->vacuum_id);
            }
        }


        const int surface_bdr_quads[6][4] = {
            {12, 13, 15, 14},
            {13, 9, 11, 15},
            {9, 8, 10, 11},
            {8, 12, 14, 10},
            {8, 9, 13, 12},
            {14, 15, 11, 10}
        };

        for (const auto& bdr: surface_bdr_quads) {
            mesh->AddBdrQuad(bdr, config->surface_bdr_id);
        }

        if (config->include_external_domain) {
            const int inf_bdr_quads[6][4] = {
                {16, 17, 21, 20},
                {17, 19, 23, 21},
                {19, 18, 22, 23},
                {18, 16, 20, 22},
                {18, 19, 17, 16},
                {20, 21, 23, 22}
            };

            for (const auto& bdr: inf_bdr_quads) {
                mesh->AddBdrQuad(bdr, config->inf_bdr_id);
            }
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
