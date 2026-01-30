#include "stroid/utils/mesh_utils.h"
#include "mfem.hpp"
#include <print>

namespace stroid::utils {
    void MarkFlippedElements(mfem::Mesh& mesh) {
        size_t total_flipped = 0;
        size_t total_elements = mesh.GetNE();
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
                total_flipped++;
            }


        }
        std::println("Marked {}/{} elements as flipped.", total_flipped, total_elements);
    }

    void MarkFlippedBoundaryElements(mfem::Mesh& mesh) {
        size_t total_flipped = 0;
        size_t total_boundary_elements = mesh.GetNBE();
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
                total_flipped++;
            }
        }
        std::println("Marked {}/{} boundary elements as flipped.", total_flipped, total_boundary_elements);
    }
}