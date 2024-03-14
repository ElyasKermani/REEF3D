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

#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"
 
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include <iostream>
#include <fstream>

using namespace chrono;
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
    floorBody->SetPos(ChVector<>(0, -1, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(floorBody);

    //Shared contact mat for all meshes

    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    // Load the .stl mesh file for the rigid body
    // Create the rigid body based on the .stl mesh

    /*
    auto body = chrono_types::make_shared<ChBody>();
    sys.Add(body);
    body->SetMass(10);
    body->SetPos(ChVector<>(0,1,0));
    body->SetInertiaXX(ChVector<>(1,1,1));
    */

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("C:/workspace/chrono_build/bin/data/models/reef3d/solid_bin.stl");

    //trimesh->Transform(ChVector<>(0,0,0), ChMatrix33<>(1.2));
    //trimesh->RepairDuplicateVertexes(1e-9);

    /*
    // compute mass inertia from mesh   REEF3D-6DOF-0-000000.stl
    double mass;
    ChVector<> cog;
    ChMatrix33<> inertia;
    double density = 1000;
    trimesh->ComputeMassProperties(true, mass, cog, inertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector<> principal_I;
    ChInertiaUtils::PrincipalInertia(inertia, principal_I, principal_inertia_rot);
    */
    
    auto vish_m = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    vish_m->SetMesh(trimesh);
    vish_m->SetMutable(false);
    vish_m->SetColor(ChColor(1.0f, 0.5f, 0.5f));
    vish_m->SetBackfaceCull(true);
    auto vis_model = chrono_types::make_shared<ChVisualModel>();
    vis_model->AddShape(vish_m);
    //body->AddVisualShape(mesh_a);

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);
    
    /*
    for (int j = 0; j < 5; ++j) {
        auto falling = chrono_types::make_shared<ChBodyAuxRef>();

        falling->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

        falling->SetMass(mass * density);
        falling->SetInertiaXX(density * principal_I);
        
        falling->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(-0.9 + ChRandom() * 1.4, 0.4 + j * 0.12, -0.9 + ChRandom() * 1.4)));

        sys.Add(falling);

        falling->AddVisualModel(vis_model);
        falling->AddCollisionShape(colli_shape);
        falling->SetCollide(true);
    }
    */

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    //falling->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

    // Set inertia
    //falling->SetMass(mass * density);
    //falling->SetInertiaXX(density * principal_I);

    // Set the absolute position of the body:
    //falling->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(0,2,0)));


    auto falling = chrono_types::make_shared<ChBody>();
    falling->SetMass(10);
    //falling->SetInteriaXX(ChVector<>(20,20,20));
    falling->SetPos(ChVector<>(0,10,0));
    
    sys.Add(falling);

    falling->AddVisualModel(vis_model);
    falling->AddCollisionShape(colli_shape);
    falling->SetCollide(true);


    /*
   // Set the triangle mesh as a collision model
    body->GetCollisionModel()->ChCollisionModel::Clear();
    body->GetCollisionModel()->ChCollisionModel::AddTriangleMesh(surfacemat, trimesh, false, false, VNULL, ChMatrix33<>(1), 0.01);
    body->GetCollisionModel()->ChCollisionModel::BuildModel();
    body->SetCollide(true);
    */

   auto cont_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
   cont_mat->SetFriction(0.2f);



   /*//Falling rigid body
    auto body = chrono_types::make_shared<ChBodyEasyBox>(7, 7.5, 0.2,      // fig size
                                                            1000,      // density
                                                            true,      // visualization?
                                                            true,      // collision?
                                                            cont_mat);  // contact material
    body->SetPos(ChVector<>(-0.5 + ChRandom() * 1, 1.4, -0.5 + ChRandom()));
    sys.Add(body);
    body->AddVisualShape(vish_m);
    */

   for (int bi = 0; bi < 20; bi++) {
        auto sphereBody = chrono_types::make_shared<ChBodyEasySphere>(0.05,      // radius size
                                                                      1000,      // density
                                                                      true,      // visualization?
                                                                      true,      // collision?
                                                                      cont_mat);  // contact material
        sphereBody->SetPos(ChVector<>(-0.5 + ChRandom() * 1, 1.4, -0.5 + ChRandom()));
        sphereBody->GetVisualShape(0)->SetColor(ChColor(0.3f, 0.3f, 0.6f));
        sys.Add(sphereBody);
    }

    /*
    // Visualization with ChVisualShapeFEA
    auto visualizemeshA = chrono_types::make_shared<ChVisualShapeFEA>(vish_m);
    visualizemeshA->SetColorscaleMinMax(0.0, 5.50);
    visualizemeshA->SetShrinkElements(true, 0.85);
    visualizemeshA->SetSmoothFaces(true);
    vish_m->AddVisualShapeFEA(visualizemeshA);
    */

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

//differnce between stl and obj
//wavefront obj?