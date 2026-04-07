#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"

#include <iostream>

namespace stroid::topology {
    void PromoteToHighOrder(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        const auto* fec = new mfem::H1_FECollection(config->order.value(), mesh.Dimension());
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
        const int nElem = mesh.GetNE();

        std::vector<bool> processed(nDofs, false);
        mfem::Array<int> vdofs;
        mfem::Vector pos(vDim);

        for (int elemID = 0; elemID < nElem; ++elemID) {
            const int attrID = mesh.GetAttribute(elemID);
            fes->GetElementVDofs(elemID, vdofs);

            for (int dofID = 0; dofID < vdofs.Size(); ++dofID) {
                const int vDof = vdofs[dofID];
                const int scalar_dof = (fes->GetOrdering() == mfem::Ordering::byNODES) ? vDof / vDim : vDof % nDofs;

                if (processed[scalar_dof]) {
                    continue; // Skip already processed dofs. This avoids doing multiple transformations of a node if it was already transformed by a neighbor
                }

                for (int d = 0; d < vDim; ++d) {
                    pos(d) = nodes(fes->DofToVDof(scalar_dof, d));
                }

                TransformPoint(pos, config, attrID);

                for (int d = 0; d < vDim; ++d) {
                    nodes(fes->DofToVDof(scalar_dof, d)) = pos(d);
                }

                processed[scalar_dof] = true;
            }
        }

    }
}
