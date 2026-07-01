#include "stroid/utils/mesh_utils.h"
#include "mfem.hpp"
#include <print>

#include "stroid/topology/curvilinear.h"

namespace stroid::utils {
    void MarkFlippedElements(mfem::Mesh& mesh) {
        for (int i = 0; i < mesh.GetNE(); i++) {
            mfem::ElementTransformation *T = mesh.GetElementTransformation(i);

            const mfem::IntegrationRule &ir = mfem::IntRules.Get(T->GetGeometryType(), 2 * T->Order());

            bool is_flipped = false;
            for (int j = 0; j < ir.GetNPoints(); j++) {
                T->SetIntPoint(&ir.IntPoint(j));
                if (T->Jacobian().Det() < 0.0) {
                    is_flipped = true;
                    break;
                }
            }

            if (is_flipped) {
                mesh.SetAttribute(i, 999);
            }
        }
    }

    void MarkFlippedBoundaryElements(mfem::Mesh& mesh) {
        for (int i = 0; i < mesh.GetNBE(); i++) {
            mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
            const mfem::IntegrationRule &ir = mfem::IntRules.Get(T->GetGeometryType(), 2 * T->Order());

            bool is_flipped = false;
            for (int j = 0; j < ir.GetNPoints(); j++) {
                T->SetIntPoint(&ir.IntPoint(j));
                const mfem::DenseMatrix &J = T->Jacobian();
                mfem::Vector pos;
                T->Transform(ir.IntPoint(j), pos);

                const double nx = J(1,0) * J(2,1) - J(2,0) * J(1,1);
                const double ny = J(2,0) * J(0,1) - J(0,0) * J(2,1);
                const double nz = J(0,0) * J(1,1) - J(1,0) * J(0,1);

                if (nx * pos(0) + ny * pos(1) + nz * pos(2) < 0.0) {
                    is_flipped = true;
                    break;
                }
            }

            if (is_flipped) {
                mesh.SetBdrAttribute(i, 500);
            }
        }
    }

    void ExportJacobianRadialProfile(mfem::Mesh& mesh, const std::string& filename) {
        std::ofstream ofs(filename);

        if (!ofs.good()) {
            throw std::runtime_error(std::format("Stroid: Could not open file {} for writing Jacobian radial profile", filename));
        }

        ofs << "Radius,DetJ,Attribute,ElementID\n";
        ofs.precision(10);

        const int sample_order = 2 * mesh.GetNodes()->FESpace()->GetMaxElementOrder() + 2;
        for (int i = 0; i < mesh.GetNE(); ++i) {
            mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
            const int attr = mesh.GetAttribute(i);

            const mfem::IntegrationRule &ir = mfem::IntRules.Get(T->GetGeometryType(), sample_order);

            for (int j = 0; j < ir.GetNPoints(); ++j) {
                T->SetIntPoint(&ir.IntPoint(j));

                mfem::Vector pos;
                T->Transform(ir.IntPoint(j), pos);

                const double r = pos.Norml2();
                const double detJ = T->Jacobian().Det();

                ofs << r << "," << detJ << "," << attr << ',' << i << "\n";
            }
        }
        ofs.close();
        std::println("Jacobian radial profile exported to {}", filename);
    }

    std::unique_ptr<mfem::Mesh> BuildProjected(const mfem::Mesh& reference, const fourdst::config::Config<config::MeshConfig>& cfg) {
        auto projected = std::make_unique<mfem::Mesh>(reference);
        topology::PromoteToHighOrder(*projected, cfg);
        topology::ProjectMesh(*projected, cfg);
        return projected;
    }
}
