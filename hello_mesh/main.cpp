#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

static const char *programName = "hello_mesh";

void OnGui()
{
    ImGui::Text(programName);
}

int main()
{
    polyscope::options::programName = programName;
    polyscope::init();
    polyscope::state::userCallback = OnGui;
    {
        std::unique_ptr<SurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = readSurfaceMesh("../data/bunny.obj");

        auto meshView = polyscope::registerSurfaceMesh("bunny", geometry->inputVertexPositions, mesh->getFaceVertexList());
        meshView->setSmoothShade(true);

        geometry->requireVertexNormals();
        meshView->addVertexVectorQuantity("normal", geometry->vertexNormals)->setEnabled(true);
    }
    polyscope::show();

    return 0;
}
