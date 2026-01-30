#pragma once

#include "mfem.hpp"

namespace stroid::utils {
    void MarkFlippedElements(mfem::Mesh& mesh);
    void MarkFlippedBoundaryElements(mfem::Mesh& mesh);
}