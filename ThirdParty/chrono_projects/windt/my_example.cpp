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

#include "chrono/physics/ChLinkRevolute.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"

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
    
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(5, 0, 5,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>(0, -1, 0)); 
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);
    
    //Shared contact mat for all meshes

    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    // Load the .stl mesh file for the rigid body
    // Create the rigid body based on the .stl mesh

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("C:/workspace/chrono_build/bin/data/models/reef3d/housing_stl.stl");

    ////

    auto trimesh2 = ChTriangleMeshConnected::CreateFromSTLFile("C:/workspace/chrono_build/bin/data/models/reef3d/blades_stl.stl");

    ChQuaternion<> rotation1;
    rotation1.Q_from_AngAxis(-CH_C_PI / 2, ChVector<>(1, 0, 0));  // 1: rotate 90� on X axis
    ChQuaternion<> rotation2;
    rotation2.Q_from_AngAxis(2 * CH_C_PI, ChVector<>(0, 1, 0));  // 2: rotate 180� on vertical Y axis
    ChQuaternion<> tot_rotation = rotation2 % rotation1;     // rotate on 1 then on 2, using quaternion product
    ChFrameMoving<> root_frame(ChVector<>(0, 0, 0), tot_rotation);
    
    auto vish_m = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    vish_m->SetMesh(trimesh);
    vish_m->SetMutable(false);
    vish_m->SetColor(ChColor(1.0f, 0.5f, 0.5f));
    vish_m->SetBackfaceCull(true);
    auto vis_model = chrono_types::make_shared<ChVisualModel>();
    vis_model->AddShape(vish_m);

    ////
    
    auto vish_b = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    vish_b->SetMesh(trimesh2);
    vish_b->SetMutable(false);
    vish_b->SetColor(ChColor(1.0f, 0.5f, 0.5f));
    vish_b->SetBackfaceCull(true);
    auto vis_model2 = chrono_types::make_shared<ChVisualModel>();
    vis_model2->AddShape(vish_b);

    // Create a shared collision shape

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.0005);

    ////

    auto colli_shape2 = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh2, false, false, 0.0005);

    

    auto housing = chrono_types::make_shared<ChBody>();
    //housing->SetMass(10);
    //housing->SetInteriaXX(ChVector<>(20,20,20));
    housing->SetBodyFixed(true);
    housing->SetPos(ChVector<>(0,0,0));
    //housing->ConcatenatePreTransformation(root_frame);

    // Rotation of -90 degrees around the X axis to make the Z-axis vertical
    ChQuaternion<> yrotation = Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0));  

    // Apply this rotation to the housing
    housing->SetRot(yrotation);

    // Determine the angle by which you want to rotate around the x-axis
    ChQuaternion<> rotation_around_y = Q_from_AngAxis(CH_C_PI, ChVector<>(1, 0, 0));

    // Apply this rotation around the x-axis to the housing
    // Make sure this is done after the previous rotation that aligned the x-axis vertically
    housing->SetRot(housing->GetRot() % rotation_around_y);


    // Now the housing's Z-axis should be aligned with the global Y-axis

    sys.Add(housing);

    housing->AddVisualModel(vis_model);
    housing->AddCollisionShape(colli_shape);
    housing->SetCollide(true);

    ////

    //house = chrono_types::make_shared<ChBodyEasyCylinder>(geometry::ChAxis::X, 0.5, 0.1, 1000, mesh_mat);
    //house->SetPos();

    auto blade = chrono_types::make_shared<ChBody>();
    blade->SetMass(1.0);
    //blade->SetInertiaXX(ChVector<>(1.0,1.0,1.0));
    //blade->SetRot(Q_from_AngAxis(CH_C_PI / 2, {0, 1, 0}));
    //blade->SetBodyFixed(true);
    blade->SetPos(ChVector<>(-0.03,-0.083,-0.009));  //-0.03,-0.082,-0.009
    //blade->ConcatenatePreTransformation(root_frame);

    sys.Add(blade);

    blade->AddVisualModel(vis_model2);
    blade->AddCollisionShape(colli_shape2);
    blade->SetCollide(true);
   
    // Assuming the axis of rotation is the Z-axis, change as needed
    //ChVector<> axis(1, 0, 0);
    //ChVector<> pivot = housing->GetPos(); // The pivot point is the position of the housing
    //ChFrame<double> joint_frame(pivot, QUNIT);  // The joint frame is the same as the pivot point

    //ChQuaternion<> rot = QUNIT;

    //ChVector<> pivot(0, 0, 1);  // Assuming the pivot is at the origin and rotation is around the Z-axis
    //ChQuaternion<> rot = Q_from_AngAxis(CH_C_PI / 2, VECT_X);  // This is just an example, adjust the angle and axis as needed
    //ChFrame<> joint_frame(pivot, rot);

    //joint->Initialize(blade, housing, joint_frame);


   // Create the revolute joint between the housing and bladedes
    //auto joint = chrono_types::make_shared<ChLinkRevolute>();
    //joint->Initialize(housing, blade, joint_frame);
    //sys.AddLink(joint);

    //revoluteJoint->Initialize(housing, blade, ChCoordsys<>(pivot, Q_from_AngAxis(CH_C_PI / 2, axis)));
    //revoluteJoint->Initialize(housing, blade, ChFrame<double>(pivot, Q_from_AngAxis(CH_C_PI / 2, axis)));
    //sys.AddLink(revoluteJoint);

    
     // Optional: Add a motor to control the bladedes' rotation
    auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->Initialize(blade, housing, ChFrame<>(0,0,0));

    sys.Add(motor);

    motor->SetSpeedFunction(chrono_types::make_shared<ChFunction_Const>(1.0));  // Set the desired rotation speed in radians per second
    

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("STL Mesh Simulation");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(0, 0.5, 0));
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

//
