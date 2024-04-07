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

#include "chrono_postprocess/ChGnuPlot.h"

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

class MyChLoaderPressure : public ChLoaderUVdistributed {
  private:
    double pressure;
    bool is_stiff;
    int num_integration_points;
    double frequency;  // frequency of the sinusoidal wave
    double time;       // current time

  public:
    MyChLoaderPressure(std::shared_ptr<ChLoadableUV> mloadable)
        : ChLoaderUVdistributed(mloadable), is_stiff(false), num_integration_points(1), time(0.0), frequency(1.0) {}

    virtual void ComputeF(const double U,        ///< parametric coordinate in surface
                          const double V,        ///< parametric coordinate in surface
                          ChVectorDynamic<>& F,  ///< Result F vector here, size must be = n.field coords.of loadable
                          ChVectorDynamic<>* state_x,  ///< if != 0, update state (pos. part) to this, then evaluate F
                          ChVectorDynamic<>* state_w   ///< if != 0, update state (speed part) to this, then evaluate F
                          ) override {
        ChVector<> mnorm = this->loadable->ComputeNormal(U, V);
        double pressure_time = pressure * sin(2 * CH_C_PI * frequency * time);
        F.segment(0, 3) = -pressure_time * mnorm.eigen();
    }

    void SetPressure(double mpressure) { pressure = mpressure; }
    double GetPressure() { return pressure; }

    void SetFrequency(double mfreq) { frequency = mfreq; }
    double GetFrequency() { return frequency; }

    void SetTime(double mtime) { time = mtime; }  // setter for time
    double GetTime() { return time; }  // getter for time

    void SetIntegrationPoints(int val) { num_integration_points = val; }
    virtual int GetIntegrationPointsU() override { return num_integration_points; }
    virtual int GetIntegrationPointsV() override { return num_integration_points; }

    void SetStiff(bool val) { is_stiff = val; }
    virtual bool IsStiff() override { return is_stiff; }
};


int main(int argc, char* argv[]) {
    
    // Create a Chrono::Engine physical system + collision system
    ChSystemNSC sys;

    //sys.Set_G_acc(ChVector<>(0, 0, 0));

    SetChronoDataPath(CHRONO_DATA_DIR);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    /*
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
    /*
    if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(0))) {
            mnode->SetFixed(true);
    }
    */

   if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(2))) {
            mnode->SetFixed(false);
    }

    sys.Add(mesh);

    // Iterate over all nodes in the mesh
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Get the position of the node
            ChVector<> pos = node->GetPos();

            // Print the position
            std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
        }
    }

    // Check if the mesh is correctly loaded
    std::cout << "Number of nodes: " << mesh->GetNnodes() << std::endl;
    std::cout << "Number of elements: " << mesh->GetNelements() << std::endl;

    /*
    // Check if the nodes are correctly fixed
    for (int i = 0; i < 20; i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            node->SetFixed(true);

            // Print the position
            ChVector<> pos = node->GetPos();
            std::cout << "Fixed node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
        }
    }
    */

    auto loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    // Set the pressure value
    double pressureValue = -120.0;  // replace with your actual pressure value

    // Iterate over all elements in the mesh
    for (int i = 0; i < mesh->GetNelements(); i++) {
        if (auto el = std::dynamic_pointer_cast<ChElementShellBST>(mesh->GetElement(i))) {
            // Create a ChLoaderPressure for the element
            auto pressure_load = chrono_types::make_shared<ChLoad<MyChLoaderPressure>>(el);
            pressure_load->loader.SetPressure(pressureValue);
            pressure_load->loader.SetStiff(true); 
            pressure_load->loader.SetIntegrationPoints(2);
            pressure_load->loader.SetFrequency(1.0);
            pressure_load->loader.SetTime(sys.GetChTime());
            loadcontainer->Add(pressure_load);

            // Print the pressure value
            std::cout << "Pressure value for element " << i << ": " << pressureValue << std::endl;
        }
    }
    
    sys.Add(loadcontainer);

    // Visualization of the FEM mesh.
    auto vis_shell_mesh = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    vis_shell_mesh->SetFEMdataType(ChVisualShapeFEA::DataType::SURFACE);
    vis_shell_mesh->SetWireframe(true);
    vis_shell_mesh->SetShellResolution(2);
    vis_shell_mesh->SetBackfaceCull(true);
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
        for (auto& load : loadcontainer->GetLoadList()) {
        if (auto pressure_load = std::dynamic_pointer_cast<ChLoad<MyChLoaderPressure>>(load)) {
            pressure_load->loader.SetTime(sys.GetChTime());
            }
        }

        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.005);

        // Print node positions
        for (int i = 0; i < mesh->GetNnodes(); i++) {
            if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
                ChVector<> pos = node->GetPos();
                std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
            }
        }
    }

    return 0;
}

//to fix 1/4 window irrlicht issue: <key>NSHighResolutionCapable</key>
//                                  <string>No</string>