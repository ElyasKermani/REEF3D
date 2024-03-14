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
// Authors: Alessandro Tasora
// =============================================================================
//
//  Demo code about using motors to impose rotation or translation between parts
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChLinkMotorRotationTorque.h"
#include "chrono/physics/ChLinkMotorRotationDriveline.h"
#include "chrono/physics/ChLinkMotorLinearPosition.h"
#include "chrono/physics/ChLinkMotorLinearSpeed.h"
#include "chrono/physics/ChLinkMotorLinearForce.h"
#include "chrono/physics/ChLinkMotorLinearDriveline.h"
#include "chrono/physics/ChShaftsMotorSpeed.h"
#include "chrono/physics/ChShaftsMotorAngle.h"
#include "chrono/physics/ChShaftsPlanetary.h"
#include "chrono/physics/ChShaftsGear.h"

#include "chrono/core/ChRealtimeStep.h"
#include "chrono/motion_functions/ChFunction_Sine.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::irrlicht;

ChCollisionSystem::Type collision_type = ChCollisionSystem::Type::BULLET;

// Shortcut function that creates two bodies (a slider and a guide) in a given position,
// just to simplify the creation of multiple linear motors in this demo.
void CreateSliderGuide(std::shared_ptr<ChBody>& guide,
                       std::shared_ptr<ChBody>& slider,
                       std::shared_ptr<ChMaterialSurface> material,
                       ChSystem& sys,
                       const ChVector<> mpos) {
    guide = chrono_types::make_shared<ChBodyEasyBox>(4, 0.3, 0.6, 1000, material);
    guide->SetPos(mpos);
    guide->SetBodyFixed(true);
    sys.Add(guide);

    slider = chrono_types::make_shared<ChBodyEasyBox>(0.4, 0.2, 0.5, 1000, material);
    slider->SetPos(mpos + ChVector<>(0, 0.3, 0));
    slider->GetVisualShape(0)->SetColor(ChColor(0.6f, 0.6f, 0.0f));
    sys.Add(slider);

    auto obstacle = chrono_types::make_shared<ChBodyEasyBox>(0.4, 0.4, 0.4, 8000, material);
    obstacle->SetPos(mpos + ChVector<>(1.5, 0.4, 0));
    obstacle->GetVisualShape(0)->SetColor(ChColor(0.2f, 0.2f, 0.2f));
    sys.Add(obstacle);
}

// Shortcut function that creates two bodies (a stator and a rotor) in a given position,
// just to simplify the creation of multiple linear motors in this demo
// (skip this and go to main() for the tutorial)

void CreateStatorRotor(std::shared_ptr<ChBody>& stator,
                       std::shared_ptr<ChBody>& rotor,
                       std::shared_ptr<ChMaterialSurface> material,
                       ChSystem& sys,
                       const ChVector<>& mpos) {
    stator =
        chrono_types::make_shared<ChBodyEasyCylinder>(geometry::ChAxis::Y, 0.5, 0.1, 1000, material);
    stator->SetPos(mpos);
    stator->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
    stator->SetBodyFixed(true);
    sys.Add(stator);

    rotor = chrono_types::make_shared<ChBodyEasyBox>(1, 0.1, 0.1, 1000, material);
    rotor->SetPos(mpos + ChVector<>(0.5, 0, -0.15));
    rotor->GetVisualShape(0)->SetColor(ChColor(0.6f, 0.6f, 0.0f));
    sys.Add(rotor);
}

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Create a ChronoENGINE physical system
    ChSystemNSC sys;
    sys.SetCollisionSystemType(collision_type);

    // Contact material shared among all objects
    auto material = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    // Create a floor that is fixed (that is used also to represent the absolute reference)
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(20, 2, 20, 3000, material);
    floorBody->SetPos(ChVector<>(0, -2, 0));
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);

    // In the following we will create different types of motors
    // - rotational motors: examples A.1, A.2, etc.
    // - linear motors, examples B.1, B.2 etc.

    // EXAMPLE A.1
    //
    // - class:   ChLinkMotorRotationSpeed
    // - type:    rotational motor
    // - control: impose a time-dependent speed=v(t)
    //
    // This is a simple type of rotational actuator. It assumes that
    // you know the exact angular speed of the rotor respect to the stator,
    // as a function of time:   angular speed = w(t).
    // Use this to simulate fans, rotating cranks, etc.
    // Note: this is a rheonomic motor that enforces the motion
    // geometrically; no compliance is allowed, this means that if the
    // rotating body hits some hard contact, the solver might give unpredictable
    // oscillatory or diverging results because of the contradiction.

    ChVector<> positionA1(-3, 2, -3);
    std::shared_ptr<ChBody> stator1;
    std::shared_ptr<ChBody> rotor1;
    CreateStatorRotor(stator1, rotor1, material, sys, positionA1);

    // Create the motor
    auto rotmotor1 = chrono_types::make_shared<ChLinkMotorRotationSpeed>();

    // Connect the rotor and the stator and add the motor to the system:
    rotmotor1->Initialize(rotor1,                // body A (slave)
                          stator1,               // body B (master)
                          ChFrame<>(positionA1)  // motor frame, in abs. coords
    );
    sys.Add(rotmotor1);

    // Create a ChFunction to be used for the ChLinkMotorRotationSpeed
    auto mwspeed =
        chrono_types::make_shared<ChFunction_Const>(CH_C_PI_2);  // constant angular speed, in [rad/s], 1PI/s =180°/s
    // Let the motor use this motion function:
    rotmotor1->SetSpeedFunction(mwspeed);

    // The ChLinkMotorRotationSpeed contains a hidden state that performs the time integration
    // of the angular speed setpoint: such angle is then imposed to the
    // constraint at the positional level too, thus avoiding angle error
    // accumulation (angle drift). Optionally, such positional constraint
    // level can be disabled as follows:
    //
    // rotmotor1->SetAvoidAngleDrift(false);

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Motors");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(1, 3, -7));
    vis->AddTypicalLights();
    vis->AddLightWithShadow(ChVector<>(20.0, 35.0, -25.0), ChVector<>(0, 0, 0), 55, 20, 55, 35, 512,
                            ChColor(0.6f, 0.8f, 1.0f));
    vis->EnableShadows();

    // Modify some setting of the physical system for the simulation, if you want
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.SetSolverMaxIterations(50);

    double timestep = 0.005;
    ChRealtimeStepTimer realtime_timer;
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        sys.DoStepDynamics(timestep);
        realtime_timer.Spin(timestep);
    }

    return 0;
}
