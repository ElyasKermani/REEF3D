// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// A very simple example that can be used as template project for
// a Chrono::Engine simulator with 3D view.
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChInertiaUtils.h"

/*
#include "chrono/physics/ChLoadBodyMesh.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoaderUV.h"
*/

/*
#include "chrono/fea/ChLoadContactSurfaceMesh.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/assets/ChVisualShapeFEA.h"
*/

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

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(20, 1, 20,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>(0, 0, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);

    //Shared contact mat for all meshes

    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("C:/workspace/chrono_build/bin/data/models/reef3d/solid_bin.stl");

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

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

    //

    auto wedge = chrono_types::make_shared<ChBodyAuxRef>();
    wedge->SetMass(10);
    //wedge->SetBodyFixed(true);

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    wedge->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

    // Set inertia
    wedge->SetMass(mass * density);
    wedge->SetInertiaXX(density * principal_I);

    // Set the absolute position of the body:
    wedge->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(0,0,0)));

    ChQuaternion<> rotation1 = Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0)); 

    wedge->SetRot(rotation1);
    //wedge->SetPos(ChVector<>(0,1,0));
    
    sys.Add(wedge);

    wedge->AddVisualModel(vis_model);
    wedge->AddCollisionShape(colli_shape);
    wedge->SetCollide(true);

    //



    //

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("STL Mesh Simulation");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 10, -1));
    vis->AddTypicalLights();
    

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.005);
    }

    return 0;
}
