// =============================================================================
// PROJECT CHRONO
//==============================================================================

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/timestepper/ChTimestepper.h"

#include "chrono/fea/ChElementShellBST.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChLoadContactSurfaceMesh.h"

#include "chrono/assets/ChVisualShapeFEA.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

int main(int argc, char* argv[]) {
    
    // Create a Chrono::Engine physical system + collision system
    ChSystemNSC sys;

    //sys.Set_G_acc(ChVector<>(0, 0, 0));

    /*

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
    floorBody->SetPos(ChVector<>(0, 0, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    //sys.Add(floorBody);
    */

    auto mesh = chrono_types::make_shared<ChMesh>();

    double density = 100;

    // Create a material
    auto elasticity = chrono_types::make_shared<ChElasticityKirchhoffIsothropic>(500, 0.33);
    auto material = chrono_types::make_shared<ChMaterialShellKirchhoff>(elasticity);
    material->SetDensity(density);

    ChMeshFileLoader::BSTShellFromObjFile(mesh, "/Users/weizhiwang/workspace/circle.obj", material, 0.01);
    if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(0))) {
            mnode->SetFixed(false);
    }
    

    //create a load on the nodes of the shell
    auto cont = chrono_types::make_shared<ChLoadContainer>();

    /*
    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            ChVector<> load(0, -1 * (i + 1), 0); // Load increases with node index
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

   /*
    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            double time = sys.GetChTime(); // Get the current simulation time
            double amplitude = -1;
            double frequency = 1;
            double loadY = amplitude * sin(2 * M_PI * frequency * time); // Sinusoidal load
            ChVector<> load(0, loadY, 0);
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

   /*
    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            double time = sys.GetChTime(); // Get the current simulation time
            double amplitude = -0.5;
            double frequency = 1;
            double phase_shift = 0.1 * i; // Phase shift increases with node index
            double loadY = amplitude * sin(2 * M_PI * frequency * (time - phase_shift)); // Wave-like load
            ChVector<> load(0, loadY, 0);
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

   /*
    // Iterate over all nodes in the mesh
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Get the position of the node
            ChVector<> pos = node->GetPos();

            // Print the position
            std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
        }
    }

    for (int i = 0; i < mesh->GetNelements(); i++) {
        if (auto el = std::dynamic_pointer_cast<ChElementShellBST>(mesh->GetElement(i))) {
            std::cout << "Element " << i << std::endl;
        }
    }
    */

    /*
    // Constants for the buoyancy force
    double waterDensity = 1000; // Density of water in kg/m^3
    double gravity = 9.81; // Acceleration due to gravity in m/s^2

    // Iterate over all nodes in the mesh
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Get the position of the node
            ChVector<> pos = node->GetPos();

            // Calculate the buoyancy force
            double volume = node->GetVolume(); // Replace with the actual volume of the node
            double buoyancyForce = waterDensity * volume * gravity;

            // Apply the buoyancy force in the upward direction
            ChVector<> load(0, buoyancyForce, 0);
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

   /*
    for (int i = 0; i < 20; i++) {
        if (auto nnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            ChVector<> pos = nnode->GetPos();

            std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
        }
    }

    // Iterate over the first 19 nodes in the mesh
    for (int i = 0; i < 20; i++) {
        if (auto rnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            rnode->SetFixed(true);
        }
    }
    */

    /*
    double pressure = 1; //Pressure in Pa
    std::cout << "Pressure of " << pressure << "Pa" << std::endl;
    //auto element = std::dynamic_pointer_cast<ChElementShellBST>(mesh->GetElement(0));

    //auto pressure_load = chrono_types::make_shared<ChLoad<ChLoaderPressure>>(element);

    for (int i = 0; i < mesh->GetNelements(); i++) {
        if (auto pel = std::dynamic_pointer_cast<ChElementShellBST>(mesh->GetElement(i))) {

            auto pressure_load = chrono_types::make_shared<ChLoad<ChLoaderPressure>>(pel);
            pressure_load->loader.SetPressure(pressure);
            pressure_load->loader.SetStiff(false);
            cont->Add(pressure_load);
            std::cout << "Element with pressure" << pressure_load << std::endl;

            auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(30));
            std::cout << "Node 0 motion" << mnode->GetPos().y() << std::endl;
        }
    }
    */

    sys.Add(mesh);
    sys.Add(cont);

    // Define the force to be applied
    ChVector<> force(0, 0.1, 0); // Replace with the actual force

    // Iterate over all nodes in the mesh
    for (int i = 20; i < 31; i++) {
        if (auto fnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Apply the force to the node
            fnode->SetForce(force);
        }
    }
    

    // Visualization of the FEM mesh.
    auto vis_shell_mesh = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    vis_shell_mesh->SetFEMdataType(ChVisualShapeFEA::DataType::SURFACE);
    vis_shell_mesh->SetWireframe(true);
    vis_shell_mesh->SetShellResolution(2);
    ////vis_shell_mesh->SetBackfaceCull(true);
    mesh->AddVisualShapeFEA(vis_shell_mesh);

    auto vis_shell_speed = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    vis_shell_speed->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    vis_shell_speed->SetColorscaleMinMax(0.0, 5.0);
    vis_shell_speed->SetWireframe(false);
    vis_shell_speed->SetShellResolution(3);
    mesh->AddVisualShapeFEA(vis_shell_speed);

    auto vis_shell_nodes = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    vis_shell_nodes->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    vis_shell_nodes->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
    vis_shell_nodes->SetSymbolsThickness(0.006);
    mesh->AddVisualShapeFEA(vis_shell_nodes);

     // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Obj-tri Mesh Simulation");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 0.5, -1));
    vis->AddTypicalLights();

    //MINRES solver
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.005);
    }

    return 0;
}
