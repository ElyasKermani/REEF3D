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
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"

#include "chrono/solver/ChIterativeSolverLS.h"
 
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

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
    

    auto mesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    LoadGeoFileToMesh("/Users/weizhiwang/workspace/solid_geo.geo", mesh);

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

    //Adding forces to the system

    auto mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(mloadcontainer);

    auto mforce = chrono_types::make_shared<ChLoadBodyForce>(body, ChVector<>(10000, 0, 0), false, ChVector<>(10, 1, 0));

    mloadcontainer->Add(mforce);

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


    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);
    

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(0.01);
    }

    return 0;
}
