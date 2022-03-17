#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <geometrycentral/numerical/linear_algebra_types.h>
#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

static const char *programName = "tutte_embedding";

struct AppData
{
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
} data;

// compute tutte embedding.
VertexData<Vector2> ComputeTutteEmbedding(ManifoldSurfaceMesh &mesh)
{
    size_t numBoundaryLoop = mesh.nBoundaryLoops();
    assert(numBoundaryLoop == 1 && "mesh boundary loop must be 1!");

    size_t numVertex = mesh.nVertices();
    SparseMatrix<float> A(numVertex, numVertex);
    A.reserve(numVertex * 8); // pre-allocate for better performance.
    Vector<float> b_x(numVertex);
    Vector<float> b_y(numVertex);
    b_x.setZero();
    b_y.setZero();

    // setup pinned boundary vertices.
    auto boundary = mesh.boundaryLoop(0);
    float boundaryAngleStep = 2.0f * PI / boundary.degree();
    int boundaryCounter = 0;
    for (auto e : boundary.adjacentEdges()) {
        auto v = e.secondVertex();
        b_x[v.getIndex()] = cos(boundaryAngleStep * boundaryCounter);
        b_y[v.getIndex()] = sin(boundaryAngleStep * boundaryCounter);
        ++boundaryCounter;
    }

    // setup weight matrix.
    using Element = Eigen::Triplet<float>;
    std::vector<Element> weight;
    weight.reserve(numVertex * 6); // average valence is 6.
    for (auto v : mesh.vertices()) {
        int vi = (int)v.getIndex();
        if (v.isBoundary()) {
            weight.emplace_back(vi, vi, 1.0f);
        } else {
            float valence = 0.0f;
            for (auto vv : v.adjacentVertices()) {
                int vvi = (int)vv.getIndex();
                weight.emplace_back(vi, vvi, 1.0f); // uniform weight.
                valence += 1.0f;
            }
            weight.emplace_back(vi, vi, -valence);
        }
    }
    A.setFromTriplets(weight.begin(), weight.end());
    weight.clear();

    SquareSolver<float> solver(A);
    Vector<float> u = solver.solve(b_x);
    Vector<float> v = solver.solve(b_y);
    VertexData<Vector2> uv(mesh);
    for (size_t i = 0; i < uv.size(); ++i) {
        uv[i] = { u[i], v[i] };
    }

    return uv;
}

void OnGui()
{
    ImGui::Text(programName);
    ImGui::Separator();
    if (ImGui::Button("DoWork")) {
        auto uv = ComputeTutteEmbedding(*data.mesh);
        // setup view.
        polyscope::getSurfaceMesh("face")
            ->addVertexParameterizationQuantity("uv", uv)
            ->setEnabled(true);

        double meshScale = 100.0;
        for (size_t i = 0; i < data.geometry->inputVertexPositions.size(); ++i) {
            data.geometry->inputVertexPositions[i].x = meshScale * uv[i].x;
            data.geometry->inputVertexPositions[i].y = meshScale * uv[i].y;
            data.geometry->inputVertexPositions[i].z = 0;
        }
        polyscope::registerSurfaceMesh("face_flattened", data.geometry->inputVertexPositions, data.mesh->getFaceVertexList())
            ->setEdgeWidth(1.0);
    }
}

int main()
{
    polyscope::options::programName = programName;
    polyscope::init();
    polyscope::state::userCallback = OnGui;
    {
        std::tie(data.mesh, data.geometry) = readManifoldSurfaceMesh("../data/face.obj");
        auto meshView = polyscope::registerSurfaceMesh("face", data.geometry->inputVertexPositions, data.mesh->getFaceVertexList());
        meshView->setSmoothShade(true);
    }
    polyscope::show();

    return 0;
}
