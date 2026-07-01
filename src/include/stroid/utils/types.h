#pragma once

#include "mfem.hpp"

#include "stroid/config/config.h"

#include <memory>
#include <expected>
#include <string>
#include <unordered_map>
#include <variant>

namespace stroid {
    enum class MFEM_MESH_TYPE {
        SERIAL,
        PARALLEL
    };

    struct StroidMesh {
        MFEM_MESH_TYPE type;
        std::unique_ptr<mfem::Mesh> mesh;
        std::unique_ptr<mfem::Mesh> reference_mesh;
        config::MeshConfig config;
        size_t refinement_levels;

        [[nodiscard]] std::expected<mfem::Mesh*, std::string> as_mesh() const {
            if (type == MFEM_MESH_TYPE::SERIAL) {
                return mesh.get();
            }
            return std::unexpected{"Mesh is not serial. Try calling as_par_mesh()"};
        }

        [[nodiscard]] std::expected<mfem::Mesh*, std::string> ref_as_mesh() const {
            if (type == MFEM_MESH_TYPE::SERIAL) {
                return reference_mesh.get();
            }
            return std::unexpected{"Reference mesh is not serial. Try calling as_par_mesh()"};
        }

        [[nodiscard]] std::expected<std::unordered_map<std::string, std::variant<int, double, std::string, bool>>, std::string> mesh_stats(bool use_ref_mesh = false) const {
            if (type != MFEM_MESH_TYPE::SERIAL) {
                return std::unexpected{"Mesh is not serial. Mesh stats currently only supports serial meshes."};
            }

            mfem::Mesh* umesh;
            if (use_ref_mesh) {
                umesh = reference_mesh.get();
            } else {
                umesh = mesh.get();
            }
            std::unordered_map<std::string, std::variant<int, double, std::string, bool>> mesh_stats;
            mesh_stats.emplace("num_elements", umesh->GetNE());
            mesh_stats.emplace("num_vertices", umesh->GetNV());
            mesh_stats.emplace("num_edges", umesh->GetNEdges());
            mesh_stats.emplace("num_faces", umesh->GetNFaces());
            mesh_stats.emplace("num_boundary_elements", umesh->GetNBE());
            mesh_stats.emplace("max_bdr_attribute_id", umesh->bdr_attributes.Max());
            mesh_stats.emplace("min_bdr_attribute_id", umesh->bdr_attributes.Min());
            mesh_stats.emplace("max_element_attribute_id", umesh->attributes.Max());
            mesh_stats.emplace("min_element_attribute_id", umesh->attributes.Min());

            return mesh_stats;
        }
    };
}
