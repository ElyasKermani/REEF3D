// =============================================================================
//
// Chrono FEA using ANCF Shell 3423 elements with modified pressure load - modified version of the demo_FEA_shellsANCF_3423.cpp
//
// =============================================================================

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/fea/ChElementShellANCF_3423.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/assets/ChVisualShapeFEA.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"
#include "chrono/physics/ChLoaderUV.h"


using namespace chrono;
using namespace chrono::fea;
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
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    ChSystemSMC sys;
    //sys.Set_G_acc(ChVector<>(0, 0, -9.8));

    GetLog() << "-----------------------------------------------------------\n";
    GetLog() << "-----------------------------------------------------------\n";
    GetLog() << "     ANCF Shell Elements demo with implicit integration    \n";
    GetLog() << "-----------------------------------------------------------\n";

    // Create a mesh, that is a container for groups of elements and their referenced nodes.
    auto mesh = chrono_types::make_shared<ChMesh>();
    // Geometry of the plate
    double plate_lenght_x = 5;
    double plate_lenght_y = 0.001;
    double plate_lenght_z = 5;
    // Specification of the mesh
    int numDiv_x = 5;
    int numDiv_z = 5;
    int N_x = numDiv_x + 1;
    int N_z = numDiv_z + 1;
    // Number of elements in the z direction is considered as 1
    int TotalNumElements = numDiv_x * numDiv_z;
    int TotalNumNodes = N_x * N_z;
    // For uniform mesh
    double dx = plate_lenght_x / numDiv_x;
    double dy = plate_lenght_y;
    double dz = plate_lenght_z / numDiv_z;;

    // Create and add the nodes
    for (int i = 0; i < TotalNumNodes; i++) {
        // Node location
        double loc_x = (i % N_x) * dx;
        double loc_y = (i / N_x) % N_z * dz;
        double loc_z = 0;

        // Node direction
        double dir_x = 0;
        double dir_y = 0;
        double dir_z = 1;

        // Create the node
        auto node =
            chrono_types::make_shared<ChNodeFEAxyzD>(ChVector<>(loc_x, loc_y, loc_z), ChVector<>(dir_x, dir_y, dir_z));

        node->SetMass(0);

        /*
        // Fix all nodes along the axis X=0
        if (i % (numDiv_x + 1) == 0)
            node->SetFixed(true);

        if (i % (numDiv_z + 1) == 0)
            node->SetFixed(true);
        */

        // Fix all nodes along the axes X=0, X=max, Z=0, and Z=max
        if (i % (numDiv_x + 1) == 0 || i % (numDiv_x + 1) == numDiv_x || i / (numDiv_x + 1) == 0 || i / (numDiv_x + 1) == numDiv_z)
            node->SetFixed(true);

        // Add node to mesh
        mesh->AddNode(node);
    }

    // Get a handle to the tip node.
    auto nodetip = std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(TotalNumNodes - 1));

    // Create an orthotropic material.
    // All layers for all elements share the same material.
    double rho = 500;
    ChVector<> E(1.0e7, 1.0e7, 1.0e7);  // Reduced Young's modulus
    ChVector<> nu(0.3, 0.3, 0.3);
    ChVector<> G(4.0e6, 4.0e6, 4.0e6);  // Reduced shear modulus
    auto mat = chrono_types::make_shared<ChMaterialShellANCF>(rho, E, nu, G);

    // Create the elements
    for (int i = 0; i < TotalNumElements; i++) {
        // Adjacent nodes
        int node0 = (i / numDiv_x) * N_x + i % numDiv_x;
        int node1 = (i / numDiv_x) * N_x + i % numDiv_x + 1;
        int node2 = (i / numDiv_x) * N_x + i % numDiv_x + 1 + N_x;
        int node3 = (i / numDiv_x) * N_x + i % numDiv_x + N_x;

        // Create the element and set its nodes.
        auto element = chrono_types::make_shared<ChElementShellANCF_3423>();
        element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(node0)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(node1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(node2)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(node3)));

        for (int i = 0; i < mesh->GetNnodes(); i++) {
            if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(i))) {
                // Get the position of the node
                ChVector<> pos = node->GetPos();

                // Print the position
                std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
            }
        }

        // Set element dimensions
        element->SetDimensions(dx, dz);

        // Add a single layers with a fiber angle of 0 degrees.
        element->AddLayer(dy, 0 * CH_C_DEG_TO_RAD, mat);

        // Set other element properties
        element->SetAlphaDamp(0.0);  // Structural damping for this element

        // Add element to mesh
        mesh->AddElement(element);
    }

    auto loadcontainer = chrono_types::make_shared<ChLoadContainer>();

    // Set the pressure value
    double pressureValue = -500.0;  // replace with your actual pressure value

    // Iterate over all elements in the mesh
    for (int i = 0; i < mesh->GetNelements(); i++) {
        if (auto el = std::dynamic_pointer_cast<ChElementShellANCF_3423>(mesh->GetElement(i))) {
            // Create a ChLoaderPressure for the element
            auto pressure_load = chrono_types::make_shared<ChLoad<MyChLoaderPressure>>(el);
            pressure_load->loader.SetPressure(pressureValue);
            pressure_load->loader.SetStiff(false); 
            pressure_load->loader.SetIntegrationPoints(2);
            pressure_load->loader.SetFrequency(3.0);
            pressure_load->loader.SetTime(sys.GetChTime());
            loadcontainer->Add(pressure_load);
        }
    }
    
    sys.Add(loadcontainer);

    // Add the mesh to the system
    sys.Add(mesh);

    // -------------------------------------
    // Options for visualization in irrlicht
    // -------------------------------------

    auto visualizemeshA = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizemeshA->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    visualizemeshA->SetColorscaleMinMax(0.0, 5.50);
    visualizemeshA->SetShrinkElements(true, 0.85);
    visualizemeshA->SetSmoothFaces(true);
    mesh->AddVisualShapeFEA(visualizemeshA);

    auto visualizemeshB = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizemeshB->SetFEMdataType(ChVisualShapeFEA::DataType::SURFACE);
    visualizemeshB->SetWireframe(true);
    visualizemeshB->SetDrawInUndeformedReference(true);
    mesh->AddVisualShapeFEA(visualizemeshB);

    auto visualizemeshC = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizemeshC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
    visualizemeshC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visualizemeshC->SetSymbolsThickness(0.004);
    mesh->AddVisualShapeFEA(visualizemeshC);

    auto visualizemeshD = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizemeshD->SetFEMglyphType(ChVisualShapeFEA::GlyphType::ELEM_TENS_STRAIN);
    visualizemeshD->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visualizemeshD->SetSymbolsScale(1);
    visualizemeshD->SetColorscaleMinMax(-0.5, 5);
    visualizemeshD->SetZbufferHide(false);
    mesh->AddVisualShapeFEA(visualizemeshD);

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("ANCF Shells");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector<>(7, 0, 7), ChVector<>(0.0, 0.0, 0.0));
    vis->AttachSystem(&sys);

    // ----------------------------------
    // Perform a dynamic time integration
    // ----------------------------------

    // Set up solver
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);

    /*
    // Set up integrator
    auto stepper = chrono_types::make_shared<ChTimestepperHHT>(&sys);
    sys.SetTimestepper(stepper);
    // Alternative way of changing the integrator:
    ////sys.SetTimestepperType(ChTimestepper::Type::HHT);
    ////auto stepper = std::static_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper());

    stepper->SetAlpha(-0.2);
    stepper->SetMaxiters(5);
    stepper->SetAbsTolerances(1e-2);
    stepper->SetStepControl(true);
    stepper->SetMinStepSize(1e-4);
    ////stepper->SetVerbose(true);
    */

    // Simulation loop

    while (vis->Run()) {
        // Update time for each pressure load
        for (auto& load : loadcontainer->GetLoadList()) {
            if (auto pressure_load = std::dynamic_pointer_cast<ChLoad<MyChLoaderPressure>>(load)) {
                pressure_load->loader.SetTime(sys.GetChTime());
            }
        }

        // Print node positions
        for (int i = 0; i < mesh->GetNnodes(); i++) {
            if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(mesh->GetNode(i))) {
                ChVector<> pos = node->GetPos();
                std::cout << "Node " << i << ": " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
            }
        }

        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.001);
    }

    return 0;
}