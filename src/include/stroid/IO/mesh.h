#pragma once
#include <string>
#include "mfem.hpp"

namespace stroid::IO {
    enum class VISUALIZATION_MODE : uint8_t {
        NONE,
        ELEMENT_ID,
        BOUNDARY_ELEMENT_ID
    };

    void SaveMesh(const mfem::Mesh& mesh, const std::string& filename);
    void SaveVTU(mfem::Mesh& mesh, const std::string& exportName);
    void ViewMesh(mfem::Mesh &mesh, const std::string& title, VISUALIZATION_MODE mode);
    void VisualizeFaceValence(mfem::Mesh& mesh);
}