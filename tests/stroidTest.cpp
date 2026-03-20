#include <gtest/gtest.h>

#include "fourdst/config/config.h"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/mapping.h"
#include "stroid/topology/topology.h"
#include "stroid/utils/mesh_utils.h"

#include <cmath>
#include <filesystem>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <array>
#include <limits>

namespace {

constexpr double kPi = 3.14159265358979323846;

using Config = fourdst::config::Config<stroid::config::MeshConfig>;


std::filesystem::path GetSourceRoot() {
    if (const char* env = std::getenv("MESON_SOURCE_ROOT")) {
        return {env};
    }
    return std::filesystem::current_path();
}

Config LoadConfigFromRepo(const std::filesystem::path& relative_path) {
    Config cfg;
    cfg.load((GetSourceRoot() / relative_path).string());
    return cfg;
}


bool IsFiniteMeshNodes(const mfem::Mesh& mesh) {
    const mfem::GridFunction* nodes = mesh.GetNodes();
    if (!nodes) {
        return false;
    }
    const int vdim = nodes->FESpace()->GetVDim();
    const int ndofs = nodes->FESpace()->GetNDofs();
    for (int i = 0; i < ndofs; ++i) {
        for (int d = 0; d < vdim; ++d) {
            const double val = (*nodes)(nodes->FESpace()->DofToVDof(i, d));
            if (!std::isfinite(val)) {
                return false;
            }
        }
    }
    return true;
}

std::map<int, int> CountVolumeAttributes(const mfem::Mesh& mesh) {
    std::map<int, int> counts;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        counts[mesh.GetAttribute(i)]++;
    }
    return counts;
}

std::map<int, int> CountBoundaryAttributes(const mfem::Mesh& mesh) {
    std::map<int, int> counts;
    for (int i = 0; i < mesh.GetNBE(); ++i) {
        counts[mesh.GetBdrAttribute(i)]++;
    }
    return counts;
}

mfem::Vector TransformCopy(const mfem::Vector& in, const Config& cfg, int attribute_id = 0) {
    mfem::Vector out = in;
    stroid::topology::TransformPoint(out, cfg, attribute_id);
    return out;
}

double IntegrateElementVolume(const mfem::Mesh& mesh, int element_id) {
    mfem::ElementTransformation* T = const_cast<mfem::Mesh&>(mesh).GetElementTransformation(element_id);
    const int order = std::max(2, 2 * T->Order() + 2);
    const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), order);

    double volume = 0.0;
    for (int j = 0; j < ir.GetNPoints(); ++j) {
        const mfem::IntegrationPoint& ip = ir.IntPoint(j);
        T->SetIntPoint(&ip);
        volume += ip.weight * std::abs(T->Weight());
    }
    return volume;
}

double ComputeMeshVolume(const mfem::Mesh& mesh) {
    double volume = 0.0;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        volume += IntegrateElementVolume(mesh, i);
    }
    return volume;
}

double ComputeMeshVolumeForAttributes(const mfem::Mesh& mesh, const std::set<int>& attributes) {
    double volume = 0.0;
    for (int i = 0; i < mesh.GetNE(); ++i) {
        if (!attributes.contains(mesh.GetAttribute(i))) {
            continue;
        }
        volume += IntegrateElementVolume(mesh, i);
    }
    return volume;
}

std::unique_ptr<mfem::Mesh> BuildProjectedMesh(const Config& cfg) {
    std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);
    return mesh;
}

double ComputeStellarVolumeWithDomainLFIntegrator(mfem::Mesh& mesh, const Config& cfg) {
    const int mesh_max_attr = mesh.attributes.Size() > 0 ? mesh.attributes.Max() : 0;
    const int cfg_max_attr = static_cast<int>(std::max({cfg->core_id, cfg->envelope_id, cfg->vacuum_id}));
    const int coeff_size = std::max(1, std::max(mesh_max_attr, cfg_max_attr));

    mfem::Vector attr_coeff(coeff_size);
    attr_coeff = 0.0;
    attr_coeff(static_cast<int>(cfg->core_id) - 1) = 1.0;
    attr_coeff(static_cast<int>(cfg->envelope_id) - 1) = 1.0;

    mfem::PWConstCoefficient stellar_coeff(attr_coeff);
    mfem::L2_FECollection fec(0, mesh.Dimension());
    mfem::FiniteElementSpace fes(&mesh, &fec);
    mfem::LinearForm lf(&fes);
    lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(stellar_coeff));
    lf.Assemble();
    return lf.Sum();
}

double ColumnNorm2(const mfem::DenseMatrix& J, int col) {
    double n2 = 0.0;
    for (int r = 0; r < J.Height(); ++r) {
        n2 += J(r, col) * J(r, col);
    }
    return std::sqrt(n2);
}

double ComputeHexEdgeRatio(const mfem::Mesh& mesh, int elem_id) {
    static constexpr std::array<std::array<int, 2>, 12> edges = {{
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}},
        {{4, 5}}, {{5, 6}}, {{6, 7}}, {{7, 4}},
        {{0, 4}}, {{1, 5}}, {{2, 6}}, {{3, 7}}
    }};

    const mfem::Element* e = mesh.GetElement(elem_id);
    if (e->GetNVertices() != 8) {
        return 1.0;
    }

    const int* v = e->GetVertices();
    double min_edge = std::numeric_limits<double>::infinity();
    double max_edge = 0.0;

    for (const auto& [a, b] : edges) {
        const double* pa = mesh.GetVertex(v[a]);
        const double* pb = mesh.GetVertex(v[b]);
        const double dx = pa[0] - pb[0];
        const double dy = pa[1] - pb[1];
        const double dz = pa[2] - pb[2];
        const double len = std::sqrt(dx * dx + dy * dy + dz * dz);
        min_edge = std::min(min_edge, len);
        max_edge = std::max(max_edge, len);
    }

    if (min_edge <= 0.0 || !std::isfinite(min_edge)) {
        return std::numeric_limits<double>::infinity();
    }
    return max_edge / min_edge;
}

struct ConditioningStats {
    double min_det = std::numeric_limits<double>::infinity();
    double max_det = 0.0;
    double min_scaled_jac = std::numeric_limits<double>::infinity();
    double max_stretch_ratio = 0.0;
    double max_edge_ratio = 0.0;
    int samples = 0;
};

ConditioningStats CollectConditioningStats(const mfem::Mesh& mesh, const std::set<int>& attrs) {
    ConditioningStats stats;

    for (int i = 0; i < mesh.GetNE(); ++i) {
        if (!attrs.empty() && !attrs.contains(mesh.GetAttribute(i))) {
            continue;
        }

        stats.max_edge_ratio = std::max(stats.max_edge_ratio, ComputeHexEdgeRatio(mesh, i));

        mfem::ElementTransformation* T = const_cast<mfem::Mesh&>(mesh).GetElementTransformation(i);
        const int order = std::max(2, 2 * T->Order() + 2);
        const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), order);

        for (int j = 0; j < ir.GetNPoints(); ++j) {
            const mfem::IntegrationPoint& ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);

            const mfem::DenseMatrix& J = T->Jacobian();
            const double det = T->Weight();
            const double abs_det = std::abs(det);

            const double c0 = ColumnNorm2(J, 0);
            const double c1 = ColumnNorm2(J, 1);
            const double c2 = ColumnNorm2(J, 2);

            const double denom = c0 * c1 * c2;
            const double scaled_jac = (denom > 0.0) ? (abs_det / denom) : 0.0;

            const double cmax = std::max({c0, c1, c2});
            const double cmin = std::max(1e-16, std::min({c0, c1, c2}));
            const double stretch_ratio = cmax / cmin;

            stats.min_det = std::min(stats.min_det, det);
            stats.max_det = std::max(stats.max_det, abs_det);
            stats.min_scaled_jac = std::min(stats.min_scaled_jac, scaled_jac);
            stats.max_stretch_ratio = std::max(stats.max_stretch_ratio, stretch_ratio);
            stats.samples++;
        }
    }

    return stats;
}

} // namespace

/**
 * @brief Test suite for the Stroid library
 */
class stroidTest : public ::testing::Test {};

/**
 * @brief Verifies the baseline block topology cardinalities in the no-vacuum case.
 * @details
 * Rationale: this is the fastest canary for accidental edits in block construction order,
 * vertex indexing, or boundary-face assembly.
 * Method: build the default skeleton and assert exact counts (3D, 16 vertices, 7 hexes, 6 bdr quads).
 * If this fails: inspect `stroid::topology::BuildSkeleton` in `src/lib/topology/topology.cpp`,
 * especially `add_box`, `stellar_shells`, and `surface_bdr_quads`, plus ID defaults in
 * `src/include/stroid/config/config.h`.
 */
TEST_F(stroidTest, BuildSkeleton_DefaultCounts) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
    EXPECT_EQ(mesh->GetNV(), 24);
    EXPECT_EQ(mesh->GetNE(), 13);
    EXPECT_EQ(mesh->GetNBE(), 12);
}

/**
 * @brief Verifies topology cardinalities when the external vacuum domain is enabled.
 * @details
 * Rationale: external-domain regressions usually surface first as wrong element/boundary counts.
 * Method: load `configs/test_external_domain.toml`, build skeleton, assert exact counts (24, 13, 12).
 * If this fails: inspect vacuum block creation and infinity boundary insertion in
 * `src/lib/topology/topology.cpp` (`vacuum_shells`, `inf_bdr_quads`) and config parsing path.
 */
TEST_F(stroidTest, BuildSkeleton_ExternalDomainCounts) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain.toml");
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);
    EXPECT_EQ(mesh->Dimension(), 3);
    EXPECT_EQ(mesh->GetNV(), 24);
    EXPECT_EQ(mesh->GetNE(), 13);
    EXPECT_EQ(mesh->GetNBE(), 12);
}

/**
 * @brief Confirms attribute bookkeeping for core/envelope/vacuum and surface/infinity boundaries.
 * @details
 * Rationale: physics coupling depends on stable material and boundary IDs, not just geometry.
 * Method: count attributes immediately after skeleton build and assert expected multiplicities.
 * If this fails: inspect element insertion attribute arguments in `BuildSkeleton` and verify
 * `core_id`, `envelope_id`, `vacuum_id`, `surface_bdr_id`, `inf_bdr_id` in config fixtures.
 */
TEST_F(stroidTest, BuildSkeleton_ExternalDomainAttributes) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain.toml");
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    ASSERT_NE(mesh, nullptr);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->core_id)), 1);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->envelope_id)), 6);
    EXPECT_EQ(volume_attr_counts.at(static_cast<int>(cfg->vacuum_id)), 6);

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    EXPECT_EQ(boundary_attr_counts.at(static_cast<int>(cfg->surface_bdr_id)), 6);
    EXPECT_EQ(boundary_attr_counts.at(static_cast<int>(cfg->inf_bdr_id)), 6);
}


/**
 * @brief Ensures `Finalize` performs refinement and preserves conforming topology.
 * @details
 * Rationale: `Finalize` is the topology gate before high-order projection; nonconforming output here
 * contaminates every downstream stage.
 * Method: compare element count pre/post finalize and assert `mesh.Conforming()`.
 * If this fails: inspect `stroid::topology::Finalize` in `src/lib/topology/topology.cpp`
 * (`FinalizeTopology`, orientation checks, `UniformRefinement` loop).
 */
TEST_F(stroidTest, Finalize_RefinementIncreasesElements) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    EXPECT_GT(mesh->GetNE(), initial_elements);
    EXPECT_TRUE(mesh->Conforming());
}

/**
 * @brief Checks exact 3D uniform-refinement scaling for default topology.
 * @details
 * Rationale: each hexahedron should split into 8; this catches subtle refine-loop regressions.
 * Method: run with fixed `refinement_levels=2` config and assert `NE_final = NE_initial * 8^2`.
 * If this fails: inspect refine-loop count and any topology-side early exits in `Finalize`.
 */
TEST_F(stroidTest, Finalize_DefaultRefinementScalesHexCountByEightPowerL) {
    const Config cfg = LoadConfigFromRepo("configs/test_refinement_l2.toml");

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    const int expected = initial_elements * 8 * 8;
    EXPECT_EQ(mesh->GetNE(), expected);
}

/**
 * @brief Checks exact 3D uniform-refinement scaling for external-domain topology.
 * @details
 * Rationale: refinement behavior must be independent of whether vacuum blocks are present.
 * Method: use fixed `refinement_levels=1` external config and assert `NE_final = NE_initial * 8`.
 * If this fails: inspect `Finalize` and verify external-domain elements are not excluded from refinement.
 */
TEST_F(stroidTest, Finalize_ExternalDomainRefinementScalesHexCountByEightPowerL) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain_refinement_l1.toml");

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    const int expected = initial_elements * 8;
    EXPECT_EQ(mesh->GetNE(), expected);
}

/**
 * @brief Validates conformance + attribute presence after refining external-domain meshes.
 * @details
 * Rationale: refinement must not silently drop regions or boundaries in multi-material meshes.
 * Method: finalize external mesh, check conforming status, growth in `NE`, and nonzero counts for expected IDs.
 * If this fails: inspect `Finalize` orientation/refinement calls and any attribute mutation side effects.
 */
TEST_F(stroidTest, Finalize_ExternalDomainConformingAndRefined) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain.toml");
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    const int initial_elements = mesh->GetNE();

    stroid::topology::Finalize(*mesh, cfg);

    EXPECT_TRUE(mesh->Conforming());
    EXPECT_GT(mesh->GetNE(), initial_elements);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->core_id)), 0);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->envelope_id)), 0);
    EXPECT_GT(volume_attr_counts.at(static_cast<int>(cfg->vacuum_id)), 0);

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    EXPECT_GT(boundary_attr_counts.at(static_cast<int>(cfg->surface_bdr_id)), 0);
    EXPECT_GT(boundary_attr_counts.at(static_cast<int>(cfg->inf_bdr_id)), 0);
}

/**
 * @brief Enforces strict attribute set invariants after finalize.
 * @details
 * Rationale: presence checks alone can miss rogue IDs introduced by buggy attribute rewrites.
 * Method: assert post-finalize volume/boundary attribute keys exactly match expected sets.
 * If this fails: inspect all calls to `SetAttribute` / `SetBdrAttribute` in topology + utility code,
 * notably `src/lib/topology/topology.cpp` and `src/lib/utils/mesh_utils.cpp`.
 */
TEST_F(stroidTest, Finalize_ExternalDomainKeepsOnlyExpectedMaterialAndBoundaryIDs) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain.toml");
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);

    stroid::topology::Finalize(*mesh, cfg);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    const std::set<int> expected_volume_ids = {
        static_cast<int>(cfg->core_id),
        static_cast<int>(cfg->envelope_id),
        static_cast<int>(cfg->vacuum_id)
    };
    for (const auto& [attr, count] : volume_attr_counts) {
        EXPECT_TRUE(expected_volume_ids.contains(attr));
        EXPECT_GT(count, 0);
    }
    EXPECT_EQ(volume_attr_counts.size(), expected_volume_ids.size());

    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);
    const std::set<int> expected_boundary_ids = {
        static_cast<int>(cfg->surface_bdr_id),
        static_cast<int>(cfg->inf_bdr_id)
    };
    for (const auto& [attr, count] : boundary_attr_counts) {
        EXPECT_TRUE(expected_boundary_ids.contains(attr));
        EXPECT_GT(count, 0);
    }
    EXPECT_EQ(boundary_attr_counts.size(), expected_boundary_ids.size());
}

/**
 * @brief Verifies high-order promotion actually attaches nodal data.
 * @details
 * Rationale: projection and most quality metrics are node-based; missing nodes means pipeline misuse.
 * Method: finalize then promote, assert nodes exist and are finite.
 * If this fails: inspect `stroid::topology::PromoteToHighOrder` in
 * `src/lib/topology/curvilinear.cpp` (`H1_FECollection`, `SetNodalFESpace`).
 */
TEST_F(stroidTest, PromoteToHighOrder_SetsNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);

    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Confirms projection does not produce NaN/Inf nodal coordinates.
 * @details
 * Rationale: finite nodes are the minimum numerical sanity condition for any downstream solver.
 * Method: run full pre-projection pipeline, project once, assert all node components are finite.
 * If this fails: inspect `ProjectMesh` and mapping functions in
 * `src/lib/topology/curvilinear.cpp` and `src/lib/topology/mapping.cpp`.
 */
TEST_F(stroidTest, ProjectMesh_ProducesFiniteNodes) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);

    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Unit-checks equiangular mapping against closed-form tangent expression.
 * @details
 * Rationale: this catches formula or branch edits before they propagate into mesh-wide projection.
 * Method: transform a known point and compare against explicit `tan(pi/4 * ratio)` expectations.
 * If this fails: inspect `ApplyEquiangular` in `src/lib/topology/mapping.cpp`, especially dominant-axis logic.
 */
TEST_F(stroidTest, ApplyEquiangular_BasicTransform) {
    mfem::Vector pos(3);
    pos(0) = 1.0;
    pos(1) = 0.5;
    pos(2) = -0.25;

    stroid::topology::ApplyEquiangular(pos);

    const double expected_y = 1.0 * std::tan(kPi / 4.0 * (0.5 / 1.0));
    const double expected_z = 1.0 * std::tan(kPi / 4.0 * (-0.25 / 1.0));
    EXPECT_NEAR(pos(1), expected_y, 1e-12);
    EXPECT_NEAR(pos(2), expected_z, 1e-12);
}

/**
 * @brief Verifies spheroidal flattening scales only the z component as configured.
 * @details
 * Rationale: flattening is intentionally simple and should remain easy to reason about.
 * Method: load flattening config, apply transform to z-axis point, check exact expected z.
 * If this fails: inspect `ApplySpheroidal` in `src/lib/topology/mapping.cpp` and fixture values in
 * `configs/test_flattening.toml`.
 */
TEST_F(stroidTest, ApplySpheroidal_FlattensZ) {
    const Config cfg = LoadConfigFromRepo("configs/test_flattening.toml");

    mfem::Vector pos(3);
    pos(0) = 0.0;
    pos(1) = 0.0;
    pos(2) = 10.0;

    stroid::topology::ApplySpheroidal(pos, cfg);

    EXPECT_NEAR(pos(2), 8.0, 1e-12);
}

/**
 * @brief Ensures axis-aligned points inside the core stay unchanged in this mapping regime.
 * @details
 * Rationale: axis points are symmetry anchors; any drift here is usually a serious branch regression.
 * Method: map `(1,0,0)` with default config and assert identity.
 * If this fails: inspect `TransformPoint` core-zone blend path and instability guard in
 * `src/lib/topology/mapping.cpp`.
 */
TEST_F(stroidTest, TransformPoint_AxisInsideCore_NoChange) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 1.0;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::TransformPoint(pos, cfg, 0);

    EXPECT_NEAR(pos(0), 1.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

/**
 * @brief Ensures axis-aligned envelope points remain unchanged under default spherical setup.
 * @details
 * Rationale: this verifies envelope branch consistency and avoids silent radial drift.
 * Method: map `(3,0,0)` and assert identity under default parameters.
 * If this fails: inspect `TransformPoint` envelope branch and normalized-direction reconstruction logic.
 */
TEST_F(stroidTest, TransformPoint_AxisEnvelope_NoChange) {
    const Config cfg;

    mfem::Vector pos(3);
    pos(0) = 3.0;
    pos(1) = 0.0;
    pos(2) = 0.0;

    stroid::topology::TransformPoint(pos, cfg, 0);

    EXPECT_NEAR(pos(0), 3.0, 1e-12);
    EXPECT_NEAR(pos(1), 0.0, 1e-12);
    EXPECT_NEAR(pos(2), 0.0, 1e-12);
}

/**
 * @brief Checks continuity near `r_core` and `r_star` transition surfaces.
 * @details
 * Rationale: discontinuities at interfaces destabilize high-order interpolation and integration.
 * Method: evaluate points at `r*(1±eps)` and assert mapped separation stays small in L2 norm.
 * If this fails: inspect transition formulas in `TransformPoint` (core blend, envelope mapping)
 * and any recent edits to `core_steepness` handling.
 */
TEST_F(stroidTest, TransformPoint_IsContinuousAcrossCoreAndStarInterfaces) {
    const Config cfg;
    constexpr double eps = 1e-6;

    mfem::Vector dir(3);
    dir(0) = 1.0;
    dir(1) = 0.6;
    dir(2) = -0.4;

    mfem::Vector near_core_left = dir;
    near_core_left *= cfg->r_core * (1.0 - eps);
    mfem::Vector near_core_right = dir;
    near_core_right *= cfg->r_core * (1.0 + eps);

    const mfem::Vector core_left_mapped = TransformCopy(near_core_left, cfg);
    const mfem::Vector core_right_mapped = TransformCopy(near_core_right, cfg);

    mfem::Vector diff = core_left_mapped;
    diff -= core_right_mapped;
    EXPECT_LT(diff.Norml2(), 1e-3);

    mfem::Vector near_star_left = dir;
    near_star_left *= cfg->r_star * (1.0 - eps);
    mfem::Vector near_star_right = dir;
    near_star_right *= cfg->r_star * (1.0 + eps);

    const mfem::Vector star_left_mapped = TransformCopy(near_star_left, cfg);
    const mfem::Vector star_right_mapped = TransformCopy(near_star_right, cfg);

    diff = star_left_mapped;
    diff -= star_right_mapped;
    EXPECT_LT(diff.Norml2(), 1e-3);
}

/**
 * @brief Verifies expected sign and axis-permutation symmetry in the spherical case.
 * @details
 * Rationale: symmetry violations usually indicate branch asymmetry bugs in mapping logic.
 * Method: compare mapped values for original, sign-flipped, and axis-swapped points.
 * If this fails: inspect dominant-axis branching in `ApplyEquiangular` and normalization flow in
 * `TransformPoint`.
 */
TEST_F(stroidTest, TransformPoint_RespectsSignAndAxisPermutationSymmetryWithoutFlattening) {
    const Config cfg;

    mfem::Vector p(3);
    p(0) = 4.0;
    p(1) = 2.0;
    p(2) = 1.0;

    mfem::Vector p_neg = p;
    p_neg *= -1.0;

    mfem::Vector p_swapped(3);
    p_swapped(0) = p(1);
    p_swapped(1) = p(0);
    p_swapped(2) = p(2);

    const mfem::Vector mapped = TransformCopy(p, cfg);
    const mfem::Vector mapped_neg = TransformCopy(p_neg, cfg);
    const mfem::Vector mapped_swapped = TransformCopy(p_swapped, cfg);

    EXPECT_NEAR(mapped_neg(0), -mapped(0), 1e-12);
    EXPECT_NEAR(mapped_neg(1), -mapped(1), 1e-12);
    EXPECT_NEAR(mapped_neg(2), -mapped(2), 1e-12);

    EXPECT_NEAR(mapped_swapped(0), mapped(1), 1e-12);
    EXPECT_NEAR(mapped_swapped(1), mapped(0), 1e-12);
    EXPECT_NEAR(mapped_swapped(2), mapped(2), 1e-12);
}

/**
 * @brief Smoke-tests mesh serialization to MFEM format.
 * @details
 * Rationale: I/O regressions are easy to miss during geometry-focused development.
 * Method: write a finalized mesh to temp storage and assert file exists and is non-empty.
 * If this fails: inspect `stroid::IO::SaveMesh` in `src/lib/IO/mesh.cpp` and local filesystem perms.
 */
TEST_F(stroidTest, SaveMesh_WritesFile) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);

    const std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    const std::filesystem::path mesh_path = tmp_dir / "stroid_test_mesh.mesh";

    stroid::IO::SaveMesh(*mesh, mesh_path.string());

    ASSERT_TRUE(std::filesystem::exists(mesh_path));
    EXPECT_GT(std::filesystem::file_size(mesh_path), 0u);

    std::error_code ec;
    std::filesystem::remove(mesh_path, ec);
}

/**
 * @brief End-to-end baseline pipeline test for the default domain.
 * @details
 * Rationale: validates the canonical operation order used by both library examples and CLI.
 * Method: execute full pipeline and assert non-empty, nodal, finite output mesh.
 * If this fails: check call-order assumptions and recent edits in
 * `src/lib/topology/topology.cpp` / `src/lib/topology/curvilinear.cpp`.
 */
TEST_F(stroidTest, EndToEnd_BuildFinalizePromoteProject) {
    const Config cfg;
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_GT(mesh->GetNE(), 0);
    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief End-to-end pipeline test for external-domain meshes.
 * @details
 * Rationale: confirms the same pipeline remains valid when vacuum blocks are included.
 * Method: run full external-domain pipeline and assert conforming + finite nodal output.
 * If this fails: inspect external-domain topology assembly and projection loops over mixed attributes.
 */
TEST_F(stroidTest, EndToEnd_ExternalDomainBuildFinalizePromoteProject) {
    const Config cfg = LoadConfigFromRepo("configs/test_external_domain.toml");
    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    EXPECT_TRUE(mesh->Conforming());
    EXPECT_GT(mesh->GetNE(), 0);
    EXPECT_NE(mesh->GetNodes(), nullptr);
    EXPECT_TRUE(IsFiniteMeshNodes(*mesh));
}

/**
 * @brief Verifies stellar volume invariance with respect to adding a vacuum shell.
 * @details
 * Rationale: the vacuum region should extend domain extent, not alter stellar mass volume.
 * Method: integrate core+envelope volume in both configs and compare relative difference.
 * If this fails: inspect attribute-filtered volume helpers in this file and mapping/topology changes
 * that may leak starside nodes into vacuum geometry.
 */
TEST_F(stroidTest, Volume_StellarDomainMatchesWithAndWithoutExternalDomain) {
    const Config no_external_cfg = LoadConfigFromRepo("configs/test_volume_no_external.toml");
    const Config with_external_cfg = LoadConfigFromRepo("configs/test_volume_with_external.toml");

    const std::unique_ptr<mfem::Mesh> no_external_mesh = stroid::topology::BuildSkeleton(no_external_cfg);
    stroid::topology::Finalize(*no_external_mesh, no_external_cfg);
    stroid::topology::PromoteToHighOrder(*no_external_mesh, no_external_cfg);
    stroid::topology::ProjectMesh(*no_external_mesh, no_external_cfg);

    const std::unique_ptr<mfem::Mesh> with_external_mesh = stroid::topology::BuildSkeleton(with_external_cfg);
    stroid::topology::Finalize(*with_external_mesh, with_external_cfg);
    stroid::topology::PromoteToHighOrder(*with_external_mesh, with_external_cfg);
    stroid::topology::ProjectMesh(*with_external_mesh, with_external_cfg);

    const std::set<int> stellar_attrs_no_external = {
        static_cast<int>(no_external_cfg->core_id),
        static_cast<int>(no_external_cfg->envelope_id)
    };
    const std::set<int> stellar_attrs_with_external = {
        static_cast<int>(with_external_cfg->core_id),
        static_cast<int>(with_external_cfg->envelope_id)
    };

    const double stellar_volume_no_external = ComputeMeshVolumeForAttributes(*no_external_mesh, stellar_attrs_no_external);
    const double stellar_volume_with_external = ComputeMeshVolumeForAttributes(*with_external_mesh, stellar_attrs_with_external);

    const double rel_diff = std::abs(stellar_volume_with_external - stellar_volume_no_external) /
                            std::max(stellar_volume_with_external, stellar_volume_no_external);
    EXPECT_LT(rel_diff, 5e-3);
}

/**
 * @brief Confirms total volume decomposition into stellar + vacuum components.
 * @details
 * Rationale: this explicitly checks that vacuum exclusion logic is doing real work, not a no-op.
 * Method: on external mesh, compute total, stellar-only, and vacuum-only volumes and enforce
 * additive consistency.
 * If this fails: inspect `ComputeMeshVolume*` helpers and region attribute IDs in config fixtures.
 */
TEST_F(stroidTest, Volume_ExternalMeshExcludesVacuumWhenRequested) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_with_external.toml");

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    const std::set<int> stellar_attrs = {
        static_cast<int>(cfg->core_id),
        static_cast<int>(cfg->envelope_id)
    };
    const std::set<int> vacuum_attr = {static_cast<int>(cfg->vacuum_id)};

    const double total_volume = ComputeMeshVolume(*mesh);
    const double stellar_volume = ComputeMeshVolumeForAttributes(*mesh, stellar_attrs);
    const double vacuum_volume = ComputeMeshVolumeForAttributes(*mesh, vacuum_attr);

    EXPECT_GT(vacuum_volume, 0.0);
    EXPECT_GT(total_volume, stellar_volume);
    EXPECT_NEAR(total_volume, stellar_volume + vacuum_volume, total_volume * 1e-9 + 1e-12);
}

/**
 * @brief Compares direct Jacobian-based stellar volume to analytic sphere volume.
 * @details
 * Rationale: anchors numerical integration against a closed-form reference in the spherical limit.
 * Method: integrate core+envelope using element transformations, compare to `4/3*pi*r_star^3`.
 * If this fails: inspect mapping spherical path (`flattening=0`) and quadrature order in
 * `IntegrateElementVolume`.
 */
TEST_F(stroidTest, Volume_SphericalStellarDomainMatchesAnalyticSphere) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");

    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);

    const std::set<int> stellar_attrs = {
        static_cast<int>(cfg->core_id),
        static_cast<int>(cfg->envelope_id)
    };

    const double measured_volume = ComputeMeshVolumeForAttributes(*mesh, stellar_attrs);
    const double analytic_volume = 4.0 / 3.0 * kPi * std::pow(cfg->r_star, 3.0);
    const double rel_err = std::abs(measured_volume - analytic_volume) / analytic_volume;

    EXPECT_LT(rel_err, 1e-2);
}

/**
 * @brief Repeats spherical analytic-volume check via MFEM `DomainLFIntegrator`.
 * @details
 * Rationale: independent integration machinery lowers the risk of helper-specific false confidence.
 * Method: build attribute-weighted `PWConstCoefficient` (core+envelope=1, vacuum=0), assemble
 * domain linear form, and compare against analytic sphere volume.
 * If this fails: inspect coefficient indexing (attr-1 convention), `ComputeStellarVolumeWithDomainLFIntegrator`,
 * and MFEM assembly setup in this test file.
 */
TEST_F(stroidTest, Volume_SphericalStellarDomainDomainLFIntegratorMatchesAnalyticSphere) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_spherical_with_external.toml");

    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);
    const double measured_volume = ComputeStellarVolumeWithDomainLFIntegrator(*mesh, cfg);
    const double analytic_volume = 4.0 / 3.0 * kPi * std::pow(cfg->r_star, 3.0);
    const double rel_err = std::abs(measured_volume - analytic_volume) / analytic_volume;

    EXPECT_LT(rel_err, 1e-2);
}

/**
 * @brief Enforces baseline conditioning bounds for the default projected mesh.
 * @details
 * Rationale: this guards against silent degradation in element quality that may still pass finiteness checks.
 * Method: sample Jacobian stats over quadrature points and assert positivity + distortion/stretch bounds.
 * If this fails: inspect mapping formulas in `src/lib/topology/mapping.cpp` and any changes to
 * refinement/order config used by `configs/test_volume_spherical_no_external.toml`.
 */
TEST_F(stroidTest, Conditioning_DefaultMeshHasPositiveJacobiansAndReasonableShape) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    const std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const ConditioningStats stats = CollectConditioningStats(*mesh, {});

    ASSERT_GT(stats.samples, 0);
    EXPECT_GT(stats.min_det, 1e-10);
    EXPECT_LT(stats.max_det / stats.min_det, 1e6);
    EXPECT_GT(stats.min_scaled_jac, 2e-2);
    EXPECT_LT(stats.max_stretch_ratio, 50.0);
    EXPECT_LT(stats.max_edge_ratio, 50.0);
}

/**
 * @brief Applies Jacobian conditioning checks independently to core, envelope, and vacuum regions.
 * @details
 * Rationale: global stats can hide localized failures; region-level checks make regressions diagnosable.
 * Method: collect conditioning statistics per attribute and enforce positive Jacobians + scaled-Jacobian floors.
 * If this fails: inspect region-specific mapping behavior in `TransformPoint` and verify attribute
 * assignment in `BuildSkeleton`.
 */
TEST_F(stroidTest, Conditioning_ExternalMeshPerRegionHasPositiveJacobians) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_spherical_with_external.toml");
    const std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    const ConditioningStats core_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->core_id)});
    const ConditioningStats envelope_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->envelope_id)});
    const ConditioningStats vacuum_stats = CollectConditioningStats(*mesh, {static_cast<int>(cfg->vacuum_id)});

    ASSERT_GT(core_stats.samples, 0);
    ASSERT_GT(envelope_stats.samples, 0);
    ASSERT_GT(vacuum_stats.samples, 0);

    EXPECT_GT(core_stats.min_det, 1e-10);
    EXPECT_GT(envelope_stats.min_det, 1e-10);
    EXPECT_GT(vacuum_stats.min_det, 1e-10);

    EXPECT_GT(core_stats.min_scaled_jac, 2e-2);
    EXPECT_GT(envelope_stats.min_scaled_jac, 2e-2);
    EXPECT_GT(vacuum_stats.min_scaled_jac, 1e-3);
}

/**
 * @brief Validates orientation quality via flipped-element and flipped-boundary markers.
 * @details
 * Rationale: this is a direct orientation sanity check using project utilities already used for debugging.
 * Method: run `MarkFlippedElements`/`MarkFlippedBoundaryElements` and assert sentinel attrs are absent.
 * If this fails: inspect Jacobian sign behavior and boundary normal orientation code in
 * `src/lib/utils/mesh_utils.cpp`, then trace upstream mapping changes.
 */
TEST_F(stroidTest, Conditioning_DefaultMeshHasNoFlippedElementsOrBoundaryFaces) {
    const Config cfg = LoadConfigFromRepo("configs/test_volume_spherical_no_external.toml");
    std::unique_ptr<mfem::Mesh> mesh = BuildProjectedMesh(cfg);

    stroid::utils::MarkFlippedElements(*mesh);
    stroid::utils::MarkFlippedBoundaryElements(*mesh);

    const auto volume_attr_counts = CountVolumeAttributes(*mesh);
    const auto boundary_attr_counts = CountBoundaryAttributes(*mesh);

    EXPECT_FALSE(volume_attr_counts.contains(999));
    EXPECT_FALSE(boundary_attr_counts.contains(500));
}

