
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChInertiaUtils.h"


#include "chrono/physics/ChLoadBodyMesh.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoaderUV.h"


#include "chrono/physics/ChLinkRevolute.h"
#include "chrono/fea/ChLoadContactSurfaceMesh.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/assets/ChVisualShapeFEA.h"


#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"
 
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

//#include "chrono/solver/ChIterativeSolverLS.h"

#include <iostream>
#include <fstream>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

int main(int argc, char* argv[]) {
    
    // Create a Chrono::Engine physical system + collision system
    ChSystemNSC sys;

    SetChronoDataPath(CHRONO_DATA_DIR);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    //create surface material
    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(20, 2, 20,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>(0, -2, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);
    
    //Shared contact mat for all meshes

    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto trimesh = ChTriangleMeshConnected::CreateFromWavefrontFile("/Users/weizhiwang/workspace/Chrono-test.obj");

    //trimesh->Transform(ChVector<>(0,0,0), ChMatrix33<>(1.2));
    //trimesh->RepairDuplicateVertexes(1e-9);

    
    // compute mass inertia from mesh   REEF3D-6DOF-0-000000.stl
    double mass;
    ChVector<> cog;
    ChMatrix33<> inertia;
    double density = 1000;
    trimesh->ComputeMassProperties(true, mass, cog, inertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector<> principal_I;
    ChInertiaUtils::PrincipalInertia(inertia, principal_I, principal_inertia_rot);

    //

    auto vish_m = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    vish_m->SetMesh(trimesh);
    vish_m->SetMutable(false);
    //vish_m->SetColor(ChColor(1.0f, 0.5f, 0.5f));
    vish_m->SetBackfaceCull(true);
    auto vis_model = chrono_types::make_shared<ChVisualModel>();
    vis_model->AddShape(vish_m);

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.00000001);

    //

    auto box1 = chrono_types::make_shared<ChBodyAuxRef>();
    box1->SetMass(10);
    //box1->SetBodyFixed(true);

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    box1->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

    // Set inertia
    box1->SetMass(mass * density);
    box1->SetInertiaXX(density * principal_I);

    // Set the absolute position of the body:
    //box1->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(0,0,0)));

    //ChQuaternion<> rotation1 = Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0)); 

    //box1->SetRot(rotation1);
    box1->SetPos(ChVector<>(0,0,0));
    
    sys.Add(box1);

    box1->AddVisualModel(vis_model);
    box1->AddCollisionShape(colli_shape);
    box1->SetCollide(true);
    
    /*
    auto load_container = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(load_container);

    ChVector<> force(0, -10000, 0);

    // Apply the second part of the load to box2 at its center
    auto load = chrono_types::make_shared<ChLoadBodyForce>(box3, force, true, box3->GetPos(), true);
    load_container->Add(load);
    */
    

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("REEF3D .obj Mesh Simulation");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 10, -1));
    vis->AddTypicalLights();
    

    // Simulation loop
    double amp = 5000.0;
    double freq = 1.0;
    double time = 0.0;

    /*
    while (vis->Run()) {
        // Calculate the oscillating force
        ChVector<> force(0, -amp * sin(2 * CH_C_PI * freq * time), 0);

        auto load_container = chrono_types::make_shared<ChLoadContainer>();

        sys.Add(load_container);

        // Apply the oscillating force to box3 at its center
        auto load1 = chrono_types::make_shared<ChLoadBodyForce>(box1, force, true, box1->GetPos(), true);
        load_container->Add(load1);

        auto load2 = chrono_types::make_shared<ChLoadBodyForce>(box2, force, true, box2->GetPos(), true);
        load_container->Add(load2);

        auto load3 = chrono_types::make_shared<ChLoadBodyForce>(box3, force, true, box3->GetPos(), true);
        load_container->Add(load3);

        auto load4 = chrono_types::make_shared<ChLoadBodyForce>(box4, force, true, box4->GetPos(), true);
        load_container->Add(load4);

        auto load5 = chrono_types::make_shared<ChLoadBodyForce>(box5, force, true, box5->GetPos(), true);
        load_container->Add(load5);

        auto load6 = chrono_types::make_shared<ChLoadBodyForce>(box6, force, true, box6->GetPos(), true);
        load_container->Add(load6);

        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.005);

        // Remove the load from the previous step
        sys.Remove(load_container);

        time += 0.005; // Increment the time
    }
    */

    while (vis->Run()) {
        auto load_container = chrono_types::make_shared<ChLoadContainer>();

        /*
        // Apply the oscillating force to each box at its center with a phase shift
        auto load1 = chrono_types::make_shared<ChLoadBodyForce>(box1, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + 0), 0), true, box1->GetPos(), true);
        load_container->Add(load1);

        auto load2 = chrono_types::make_shared<ChLoadBodyForce>(box2, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + CH_C_PI / 6), 0), true, box2->GetPos(), true);
        load_container->Add(load2);

        auto load3 = chrono_types::make_shared<ChLoadBodyForce>(box3, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + CH_C_PI / 3), 0), true, box3->GetPos(), true);
        load_container->Add(load3);

        auto load4 = chrono_types::make_shared<ChLoadBodyForce>(box4, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + CH_C_PI / 2), 0), true, box4->GetPos(), true);
        load_container->Add(load4);

        auto load5 = chrono_types::make_shared<ChLoadBodyForce>(box5, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + 2 * CH_C_PI / 3), 0), true, box5->GetPos(), true);
        load_container->Add(load5);

        auto load6 = chrono_types::make_shared<ChLoadBodyForce>(box6, ChVector<>(0, -amp * sin(2 * CH_C_PI * freq * time + 5 * CH_C_PI / 6), 0), true, box6->GetPos(), true);
        load_container->Add(load6);

        sys.Add(load_container);
        */

        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.005);

        sys.Remove(load_container);

        time += 0.005; // Increment the time
    }

    return 0;
}
