// =============================================================================
// PROJECT CHRONO
//==============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChInertiaUtils.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/assets/ChVisualShapeFEA.h"

#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChElementShellANCF_3423.h"
#include "chrono/fea/ChMaterialShellKirchhoff.h"


#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::geometry;
using namespace chrono::irrlicht;


void GeoFileToFEAMesh(
    std::shared_ptr<ChMesh> mesh,
    const char* filename,
    std::shared_ptr<ChMaterialShellKirchhoff> my_material,
    double my_thickness,
    ChVector<> pos_transform = ChVector<>(0, 0, 0),
    ChMatrix33<> rot_transform = ChMatrix33<>(1)
) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    std::string line;
    std::vector<std::shared_ptr<ChNodeFEAxyz>> nodes;
    std::vector<std::array<int, 3>> faces; // Assuming triangular elements

    // Skip header lines up to "NPointAttrib"
    while (std::getline(file, line) && line.find("NPointAttrib") == std::string::npos) {}

    // Read nodes
    for (int i = 0; i < 8; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        double x, y, z;
        iss >> x >> y >> z; // Skipping the last value (assumed 1 based on the provided format)
        ChVector<double> pos(x, y, z);
        pos = rot_transform * pos + pos_transform;
        auto node = chrono_types::make_shared<ChNodeFEAxyz>(pos);
        mesh->AddNode(node);
        nodes.push_back(node);
    }

    // Skip header lines
    while (std::getline(file, line) && line != "Run 12 Poly") {}

    // Read faces
    while (std::getline(file, line)) {
        if (line.empty() || line == "beginExtra") break; // Assuming faces section ends with an empty line or "beginExtra"
        std::istringstream iss(line);
        char dummy; // To skip '3' and '<'
        int i0, i1, i2;
        iss >> dummy >> dummy >> i0 >> i1 >> i2;
        faces.push_back({i0, i1, i2});
    }

    // Create elements
    for (const auto& face : faces) {
        auto element = chrono_types::make_shared<ChElementShellANCF>();
        element->SetNodes(nodes[face[0]], nodes[face[1]], nodes[face[2]]);
        element->AddLayer(my_thickness, 0 /* angle */, my_material);
        mesh->AddElement(element);
    }
}


/*
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChElementShellANCF.h"
#include "chrono/fea/ChMaterialShellANCF.h"
#include "chrono/core/ChVector.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace chrono;
using namespace chrono::fea;

// Function to read .geo file and create FEA mesh
void LoadGeoFileToFEAMesh(const std::string& filename, std::shared_ptr<ChMesh>& feaMesh, ChSystemNSC& system) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: cannot open file " + filename);
    }

    feaMesh = std::make_shared<ChMesh>();

    // Temporary storage for nodes
    std::vector<std::shared_ptr<ChNodeFEAxyz>> nodes;

    // Read nodes
    std::string line;
    while (std::getline(file, line) && line != "Run 12 Poly") {
        std::istringstream iss(line);
        double x, y, z, dummy;
        if (iss >> x >> y >> z >> dummy) {
            auto node = std::make_shared<ChNodeFEAxyz>(ChVector<>(x, y, z));
            feaMesh->AddNode(node);
            nodes.push_back(node);
        }
    }

    // Create material
    auto mat = std::make_shared<ChMaterialShellANCF>(500, 0.3, 0.01);

    // Read triangles and create shell elements
    while (std::getline(file, line)) {
        int a, b, c; // Triangle node indices
        if (sscanf(line.c_str(), "3 < %d %d %d", &a, &b, &c) == 3) {
            auto element = std::make_shared<ChElementShellANCF>();
            element->SetNodes(nodes[a], nodes[b], nodes[c], nodes[a]); // Re-use the first node for the fourth corner
            element->AddLayer(0.01, 0, mat); // Thickness, angle (radians), material
            feaMesh->AddElement(element);
        }
    }

    // Add the FEA mesh to the system
    system.Add(feaMesh);
}

int main() {
    ChSystemNSC system;

    auto feaMesh = std::make_shared<ChMesh>();

    // Example: Load .geo file and create FEA mesh
    std::string geoFilePath = "path/to/your/file.geo";
    LoadGeoFileToFEAMesh(geoFilePath, feaMesh, system);

    // Continue setting up your simulation...

    return 0;
}
*/


/*
void LoadGeoFileToMesh(const std::string& filename, std::shared_ptr<ChTriangleMeshConnected>& mesh) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    mesh = std::make_shared<ChTriangleMeshConnected>();

    // Skip header NPointAttrib
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("NPointAttrib") != std::string::npos) {
            break;
        }
    }

    // Read node definitions
    std::vector<ChVector<>> nodes;
    for (int i = 0; i < 8; ++i) { // Assumes 8 nodes as per the file structure
        std::getline(file, line);
        std::istringstream iss(line);
        double x, y, z, dummy;
        if (!(iss >> x >> y >> z >> dummy)) {
            std::cerr << "Error reading node coordinates" << std::endl;
            return;
        }
        nodes.push_back(ChVector<>(x, y, z));
    }

    // Skip to polygons definition
    while (std::getline(file, line)) {
        if (line.find("Run") != std::string::npos) {
            break;
        }
    }

    // Read polygons and create elements
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char dummy; // To skip the '3' and '<' characters
        int a, b, c; // Node indices
        if (!(iss >> dummy >> dummy >> a >> b >> c)) {
            break; // Assumes end of file or malformed line
        }
        mesh->addTriangle(nodes[a], nodes[b], nodes[c]);
    }
}
*/

int main(int argc, char* argv[]) {
    
    // Create a Chrono::Engine physical system + collision system
    ChSystemNSC sys;

    SetChronoDataPath(CHRONO_DATA_DIR);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    
    //create surface material
    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(10, 0, 10,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>(0, -1, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);

    auto mesh = chrono_types::make_shared<ChMesh>();

    auto material = chrono_types::make_shared<ChMaterialShellKirchhoff>(100, 1e8, 0.3);

    GeoFileToFEAMesh(mesh, "/Users/weizhiwang/workspace/solid_geo.geo", material, 0.01);

    sys.Add(mesh);

    // compute mass inertia from mesh
    double mass;
    ChVector<> cog;
    ChMatrix33<> inertia;
    double density = 1000;
    mesh->ComputeMassProperties(true, mass, cog, inertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector<> principal_I;
    ChInertiaUtils::PrincipalInertia(inertia, principal_I, principal_inertia_rot);
    
    //Shared contact mat for all meshes
    /*
    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, mesh, false, false, 0.0005);

    auto body = chrono_types::make_shared<ChBodyAuxRef>();
    sys.Add(body);
    //body->SetBodyFixed(true);
    //body->SetPos(ChVector<>(0, 0, 0));

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    body->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

    // Set inertia
    body->SetMass(mass * density);
    body->SetInertiaXX(density * principal_I);

    // Set the absolute position of the body:
    body->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(-5,0,0)));

    ChQuaternion<> rotation1 = Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0)); 

    body->SetRot(rotation1);

    auto mesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    mesh_shape->SetMesh(mesh);
    mesh_shape->SetMutable(false);
    mesh_shape->SetBackfaceCull(true);
    auto vis_model = chrono_types::make_shared<ChVisualModel>();
    vis_model->AddShape(mesh_shape);

    body->AddVisualModel(vis_model);
    body->AddCollisionShape(colli_shape);
    body->SetCollide(true);
    */

    //Adding forces to the system

    // -----------------------------------------------------------------
    // Set visualization of the FEM mesh.

    auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_MZ);
    mvisualizebeamA->SetColorscaleMinMax(-400, 200);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    mesh->AddVisualShapeFEA(mvisualizebeamA);

    auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    mesh->AddVisualShapeFEA(mvisualizebeamC);

    /*
    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Loads on beams");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector<>(0.5, 0.0, -3.0), ChVector<>(0.5, 0.0, 0.0));
    vis->AttachSystem(&sys);
    */

     // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("GEO Mesh Simulation");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    //vis->AddCamera(ChVector<>(0, 1, -1));
    vis->AddCamera(ChVector<>(0, 15, -20), ChVector<>(0, 0, 0));
    vis->AddTypicalLights();

    // -----------------------------------------------------------------

    // Setup a MINRES solver. For FEA one cannot use the default PSOR type solver.
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(200);
    solver->SetTolerance(1e-15);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(false);

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.01);
    }

    return 0;
}
