#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"

#include <iostream>
#include <memory>

namespace stroid::topology {
    void PromoteToHighOrder(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        const auto* fec = new mfem::H1_FECollection(config->order, mesh.Dimension());
        auto* fes = new mfem::FiniteElementSpace(&mesh, fec, mesh.SpaceDimension());
        mesh.SetNodalFESpace(fes);
    }

    void ProjectMesh(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        if (!mesh.GetNodes()) {
            std::cerr << "Error: Mesh has no nodes to project. Call PromoteToHighOrder first." << std::endl;
            return;
        }

        mfem::GridFunction& nodes = *mesh.GetNodes(); // Already confirmed not null
        const mfem::FiniteElementSpace* fes = nodes.FESpace();

        const int vDim = fes->GetVDim();
        const int nDofs = fes->GetNDofs();

        mfem::Vector pos(vDim);

        for (int i = 0; i < nDofs; ++i) {
            for (int d = 0; d < vDim; ++d) {
                pos(d) = nodes(fes->DofToVDof(i, d));
            }

            TransformPoint(pos, config, 0);

            for (int d = 0; d < vDim; ++d) {
                nodes(fes->DofToVDof(i, d)) = pos(d);
            }
        }

    }
}
