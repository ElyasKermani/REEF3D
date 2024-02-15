/*
 <DSPHCHRONOLIB>  Copyright (c) 2022 by 
 
 Iván Martínez Estévez (ivan.martinez.estevez@uvigo.es). Universidade de Vigo, Spain.
 Dr José M. Domínguez (jmdominguez@uvigo.es). Universidade de Vigo, Spain.
 Dr Ricardo Canelas (ricardo.canelas@bentley.com). Bentley Systems, Lisbon, Portugal.
 Dr Bonaventura Tagliafierro (btagliafierro@gmail.com). University of Salerno, Italy.
 Dr Orlando García Feal (orlando@uvigo.es). Universidade de Vigo, Spain.
 Professor Alejandro J.C. Crespo (alexbexe@uvigo.es). Universidade de Vigo, Spain.
 Professor Moncho Gómez Gesteira (mggesteira@uvigo.es). Universidade de Vigo, Spain. 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 
 This file is part of DSPHChronoLib. 

 DSPHChronoLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DSPHChronoLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DSPHChronoLib. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file DSPHChronoLibSC.cpp \brief Implements the class \ref DSPHChronoLibSC.

#include "DSPHChronoLib.h"
#include "Functions.h"
#include "FunDSPHChrono.h"

#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

//-------- Physics dependencies ---------
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChMaterialSurfaceNSC.h"
#include "chrono/physics/ChLinksAll.h"

//-------- Solver dependencies ---------
#include "chrono/solver/ChSolverBB.h"
#include "chrono/solver/ChSolverMINRES.h"

//-------- Core dependencies ---------
#include "chrono/core/ChVector.h"

#include <cmath>
#include <string>
#include <fstream>
#include <omp.h>
#include <memory.h>

//-Use the namespace of Chrono.
using namespace ::chrono;
using namespace ::chrono::collision;

//-------------------------------------------------------------------------------------
/// ---------------------------DSPHChronoLibSC-------------------------------------
//-------------------------------------------------------------------------------------

//==============================================================================
/// Constructor.
//==============================================================================
DSPHChronoLibSC::DSPHChronoLibSC(const JChronoData &chdata)
  : DSPHChronoLib(chdata){
  
  //-Create a Chrono physical system.
  if(UseSMC) MphysicalSystem=new ChSystemSMC;
  else       MphysicalSystem=new ChSystemNSC;

  MphysicalSystem->Set_G_acc(ChVector<>(0.0,0.0,0.0)); //-Gravity force already comming from DPSH

  if     (SolverIx==fdsch::as_integer(ChSolver::Type::BARZILAIBORWEIN)) MphysicalSystem->SetSolver(std::make_shared<ChSolverBB>());
  else if(SolverIx==fdsch::as_integer(ChSolver::Type::MINRES))          MphysicalSystem->SetSolver(std::make_shared<ChSolverMINRES>()); 
  
  MphysicalSystem->SetMaxPenetrationRecoverySpeed(0.1);
  MphysicalSystem->SetMinBounceSpeed(0.05);

  //- Cofigures Collision model
  //Setting the value of the envelope/margin (These settings will affect collision shapes that are created after these function calls.)
  collision::ChCollisionModel::SetDefaultSuggestedEnvelope(ChData.GetDp());
  collision::ChCollisionModel::SetDefaultSuggestedMargin(ChData.GetDp());
  if(UseSMC){
    //-Set default effective radius of curvature for all SCM contacts.  
    collision::ChCollisionInfo::SetDefaultEffectiveCurvatureRadius(1); //Default=0.1 
  }
}

//==============================================================================
/// Destructor.
//==============================================================================
DSPHChronoLibSC::~DSPHChronoLibSC(){
  if(MphysicalSystem)MphysicalSystem=NULL;
}

//==============================================================================
/// Loads data for bodies and configures objects.
//==============================================================================
void DSPHChronoLibSC::Config(std::string dirout,bool svdata,bool simulate2d){

  DirOut=dirout;
  Simulate2D=simulate2d;
  CollisionCoef=ChData.GetDp()*ChData.GetCollisionDp();
  
  //-Loop for chrono bodies
  for(unsigned cb=0; cb<ChData.GetBodyCount(); cb++){
    const JChBody* body=ChData.GetBody(cb);
    if(body->Type==JChBody::BD_Floating)     ConfigFloating(body);
    else if(body->Type==JChBody::BD_Fixed)   ConfigFixed(body);
    else if(body->Type==JChBody::BD_Moving)  ConfigMoving(body);
    else fun::Run_ExceptioonFun("Type of body is invalid.");
  } //all bodies are created in chrono

  //-Loop for chrono links
  for(unsigned ck=0; ck<ChData.GetLinkCount(); ck++){
    const JChLink* link=ChData.GetLink(ck);
    if(link->Type!=JChLink::LK_Hinge        &&link->Type!=JChLink::LK_Pulley
      &&link->Type!=JChLink::LK_Spheric     &&link->Type!=JChLink::LK_PointLine
      &&link->Type!=JChLink::LK_LinearSpring&&link->Type!=JChLink::LK_CoulombDamping
      &&link->Type!=JChLink::LK_PointFrame)fun::Run_ExceptioonFun("Type of link is invalid.");
    unsigned nb=link->GetBodyRefCount();
    if(nb==0)fun::Run_ExceptioonFun("At least one body per link is necessary.");
    const JChBody* body0=link->GetBodyRef(0);
    const JChBody* body1=(nb>1?link->GetBodyRef(1):NULL);

    auto fb0=MphysicalSystem->SearchBody(body0->IdName.c_str());
    auto fb1=(nb>1?MphysicalSystem->SearchBody(body1->IdName.c_str()):std::make_shared<ChBody>(UseSMC?ChMaterialSurface::SMC:ChMaterialSurface::NSC)); // - if it exists recover the body from the chrono system, if not create a new one to populate after
    if(nb<2){ //-Create a virtual hinge body, this is attached to a point in space
      fb1->SetBodyFixed(true);
      fb1->SetName("virtual_body");
      MphysicalSystem->Add(fb1);
      fb1->SetCollide(false);
    }
    //-Create the desired link
    if(link->Type==JChLink::LK_Hinge){
      const JChLinkHinge* linktype=(const JChLinkHinge*)link;
      auto MyLink=std::make_shared<ChLinkLockRevolute>();
      const tdouble3 rotpoint=linktype->GetRotPoint();
      const tdouble3 rotvector=fdsch::VecUnitary(linktype->GetRotVector());
      ChVector<> rev_point(rotpoint.x,rotpoint.y,rotpoint.z); // rotation point
      //now we need to build a quaternion that rotates the default z axis to our desired axis
      const double ang=acos(rotvector.z);// *180. / 3.14159265; //angle shortcut, since the default is always (0,0,1)
      tdouble3 to_q;
      to_q.x=rotvector.y; to_q.y=-rotvector.x; to_q.z=0; //external product shortcut, this is the normal vector to the default z from chrono and our desired rotation vector
      ChVector<> rotvec(to_q.x,to_q.y,to_q.z); // rotation vector
      rotvec.Normalize();
      MyLink->Initialize(fb0,fb1,ChCoordsys<>(rev_point,Q_from_AngAxis(ang,rotvec)));
      MyLink->SetNameString(linktype->Name);
      ChLinkForce *force=new ChLinkForce;
      force->Set_active(true);
      force->Set_K(linktype->GetStiffness()); // Torsional stiffness, to set, in [N*m/rad]
      force->Set_R(linktype->GetDamping()); // Torsional damping, to set, in [N*m*s/rad]
      MyLink->SetForce_R(force); //-Add a rotational spring damper to the revolute joint
      MphysicalSystem->AddLink(MyLink);
    }
    if(link->Type==JChLink::LK_Pulley){
      const JChLinkPulley* linktype=(const JChLinkPulley*)link;
      auto MyLink=std::make_shared<ChLinkPulley>();
      const tdouble3 rotpoint=linktype->GetRotPoint();
      const tdouble3 rotvector=fdsch::VecUnitary(linktype->GetRotVector());
      ChVector<> rev_point(rotpoint.x,rotpoint.y,rotpoint.z); // rotation point
      //now we need to build a quaternion that rotates the default z axis to our desired axis
      const double ang=acos(rotvector.z);// *180. / 3.14159265; //angle shortcut, since the default is always (0,0,1)
      tdouble3 to_q;
      to_q.x=rotvector.y; to_q.y=-rotvector.x; to_q.z=0; //external product shortcut, this is the normal vector to the default z from chrono and our desired rotation vector
      ChVector<> rotvec(to_q.x,to_q.y,to_q.z); // rotation vector
      rotvec.Normalize();
      //-Create a virtual fixed body to allow only the rotation motions
      auto fixed=std::make_shared<ChBody>(UseSMC?ChMaterialSurface::SMC:ChMaterialSurface::NSC);
      fixed->SetBodyFixed(true);
      fixed->SetName("virtual_body_pulley");
      MphysicalSystem->Add(fixed);
      fixed->SetCollide(false);
      //-Apply rotation to slave object
      //-Create revolute for slave body
      auto link_rev=std::make_shared<ChLinkLockRevolute>();
      link_rev->Initialize(fb1,fixed,ChCoordsys<>(rev_point,Q_from_AngAxis(ang,rotvec)));
      MphysicalSystem->AddLink(link_rev);
      //-Initialize the pulley between two objects
      MyLink->Initialize(fb0,fb1,ChCoordsys<>(VNULL,Q_from_AngAxis(ang,rotvec)));
      MyLink->Set_local_shaft1(ChFrame<>(VNULL,Q_from_AngAxis(ang,rotvec)));
      MyLink->Set_local_shaft2(ChFrame<>(VNULL,Q_from_AngAxis(ang,rotvec)));
      MyLink->Set_r1(linktype->GetRadius());
      MyLink->Set_r2(linktype->GetRadius2());
      MyLink->SetNameString(linktype->Name);
      MphysicalSystem->AddLink(MyLink);
    }
    if(link->Type==JChLink::LK_Spheric){
      const JChLinkSpheric* linktype=(const JChLinkSpheric*)link;
      auto MyLink=std::make_shared<ChLinkLockSpherical>();
      const tdouble3 rotpoint=linktype->GetRotPoint();
      ChVector<> rev_point(rotpoint.x,rotpoint.y,rotpoint.z); // rotation point
      MyLink->Initialize(fb0,fb1,ChCoordsys<>(rev_point));
      MyLink->SetNameString(linktype->Name);
      ChLinkForce *force=new ChLinkForce;
      force->Set_active(true);
      force->Set_K(linktype->GetStiffness()); // Torsional stiffness, to set, in [N*m/rad]
      force->Set_R(linktype->GetDamping()); // Torsional damping, to set, in [N*m*s/rad]
      MyLink->SetForce_R(force); //-Add a rotational spring damper to the revolute joint
      MphysicalSystem->AddLink(MyLink);
    }
    if(link->Type==JChLink::LK_PointLine){
      const JChLinkPointLine* linktype=(const JChLinkPointLine*)link;
      const int nlk=(linktype->GetRotVector()==TDouble3(0)?1:(linktype->GetRotVector2()==TDouble3(0)?2:3));
      for(int clk=0; clk<nlk; clk++){
        auto MyLink=std::make_shared<ChLinkLockPointLine>(); //-this stupid thing has the xx axis as default, because of reasons. 
        const tdouble3 sldvector=fdsch::VecUnitary(linktype->GetSlidingVector());
        tdouble3 rotpoint=linktype->GetRotPoint();
        if(clk==1)rotpoint=rotpoint+linktype->GetRotVector();
        if(clk==2)rotpoint=rotpoint+linktype->GetRotVector2();
        //const tdouble3 rotvector=VecUnitary(linktype->GetRotVector());
        ChVector<> rev_point(rotpoint.x,rotpoint.y,rotpoint.z); // rotation point				
        //now we need to build a quaternion that rotates the default z axis to our desired axis
        //const double ang=acos(rotvector.z);
        const double ang=acos(sldvector.x);
        tdouble3 to_q;
        to_q.x=0; to_q.y=-sldvector.z; to_q.z=-sldvector.y; //external product shortcut
        ChVector<> rotvec(to_q.x,to_q.y,to_q.z); // rotation vector
        rotvec.Normalize();
        MyLink->Initialize(fb0,fb1,ChCoordsys<>(rev_point,Q_from_AngAxis(ang,rotvec)));
        std::string linkname=linktype->Name;
        if(clk==1)linkname=linkname+"_AUTO";
        if(clk==2)linkname=linkname+"_AUTO2";
        MyLink->SetNameString(linkname);
        ChLinkForce *force=new ChLinkForce;
        force->Set_active(true);
        force->Set_K(linktype->GetStiffness()); // Torsional stiffness, to set, in [N*m/rad]
        force->Set_R(linktype->GetDamping()); // Torsional damping, to set, in [N*m*s/rad]
        MyLink->SetForce_R(force); //-Add a rotational spring damper to the lock point joint			
        MphysicalSystem->AddLink(MyLink);
      }
    }
    if(link->Type==JChLink::LK_LinearSpring){
      const JChLinkLinearSpring* linktype=(const JChLinkLinearSpring*)link;
      //MySpringForce force;
      auto MyLink=std::make_shared<ChLinkSpring>(); //-this one needes 2 points.
      const tdouble3 pointfb0=linktype->GetPointfb0();
      const tdouble3 pointfb1=linktype->GetPointfb1();
      ChVector<> point0(pointfb0.x,pointfb0.y,pointfb0.z); // point in fb0
      ChVector<> point1(pointfb1.x,pointfb1.y,pointfb1.z); // point in fb1			
      MyLink->Initialize(fb0,fb1,false,point0,point1,false,linktype->GetRestLength());
      MyLink->SetNameString(linktype->Name);
      MyLink->Set_SpringK(linktype->GetStiffness());
      MyLink->Set_SpringR(linktype->GetDamping());
      //MyLink->Set_SpringCallback(&force);
      MphysicalSystem->AddLink(MyLink);
    }
    if(link->Type==JChLink::LK_CoulombDamping){
      const JChLinkCoulombDamping* linktype=(const JChLinkCoulombDamping*)link;
      auto MyLink=std::make_shared<ChLinkCoulombDamping>(); //-this one needes 2 points.
      const tdouble3 pointfb0=linktype->GetPointfb0();
      const tdouble3 pointfb1=linktype->GetPointfb1();
      ChVector<> point0(pointfb0.x,pointfb0.y,pointfb0.z); // point in fb0
      ChVector<> point1(pointfb1.x,pointfb1.y,pointfb1.z); // point in fb1			
      MyLink->Initialize(fb0,fb1,false,point0,point1,false,linktype->GetRestLength(),linktype->GetCoulombDamping());
      MyLink->SetNameString(linktype->Name);
      MyLink->Set_SpringK(0); //Not used for CoulombDamping
      MyLink->Set_SpringR(0); //Not used for CoulombDamping
      //MyLink->Set_SpringCallback(&force);
      MphysicalSystem->AddLink(MyLink);
    }
  } //all links are created in chrono

  //-Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  if(svdata)SaveForcesHead();
  RunState=RSTATE_Init;
}

//==============================================================================
/// Configures floating bodies
//==============================================================================
void DSPHChronoLibSC::ConfigFloating(const JChBody* body){
  const double mass=body->GetMass();
  const tdouble3 center=body->GetCenter(); //not used as it is read from a .obj mesh
  ChVector<> COG(center.x,center.y,center.z);
  const tmatrix3f inertiaf=ToTMatrix3f(body->GetInertia()); //also computed from the mesh, but we want it to be coherent with DSPH inertia computation
  std::string objname=ChData.GetDataDir()+body->GetModelFile();
  //printf("Mo---------> [%s]\n",objname.c_str());
  if(!body->GetModelFile().empty()&&!fun::FileExists(objname))fun::Run_ExceptioonFileFun("Model file for floating body is missing.",objname);
  auto RigBody=std::make_shared< ChBodyEasyMesh >(objname.c_str(),1000,false,ChData.GetUseCollision(),CollisionCoef,false,(UseSMC?ChMaterialSurface::SMC:ChMaterialSurface::NSC));
  //if mesh is not read, it still creates a body object
  MphysicalSystem->Add(RigBody);  //Added to chrono scene
  ChMatrix33<> global_csys(1); //-global coord system
  RigBody->SetFrame_COG_to_REF(ChFrame<>(COG,global_csys)); //setting the COG according to DSPH, not by mesh baricenter
  //Need to attribute collision parameters				
  //Might be important to define object families to disregard contacts between mechanism parts for example
  if(ChData.GetUseCollision()) ConfigSurfaceBody(*body,RigBody.get());
  RigBody->SetMass(mass);
  RigBody->SetInertiaXX(ChVector<>(inertiaf.a11,inertiaf.a22,inertiaf.a33));
  RigBody->SetInertiaXY(ChVector<>(inertiaf.a12,inertiaf.a23,inertiaf.a13));
  RigBody->SetNameString(body->IdName);
  RigBody->SetBodyFixed(false);
  //-Apply the initial imposed velocity.
  ApplyInitialVel(*body,RigBody.get());
}

//==============================================================================
/// Configures fixed bodies
//==============================================================================
void DSPHChronoLibSC::ConfigFixed(const JChBody *body){
  //Getting fixed boundaries
  std::string objname=ChData.GetDataDir()+body->GetModelFile();
  if(!body->GetModelFile().empty()&&!fun::FileExists(objname))fun::Run_ExceptioonFileFun("Model file for fixed body is missing.",objname);
  auto RigBound=std::make_shared< ChBodyEasyMesh >(objname.c_str(),1000,false,ChData.GetUseCollision(),CollisionCoef,false,(UseSMC?ChMaterialSurface::SMC:ChMaterialSurface::NSC));
  //if mesh is not read, it still creates a body object
  MphysicalSystem->Add(RigBound);  //Added to chrono scene
  RigBound->SetNameString(body->IdName);
  RigBound->SetBodyFixed(true);
  //Need to attribute collision parameters				
  //Might be important to define object families to disregard contacts between mechanism parts for example
  if(ChData.GetUseCollision()) ConfigSurfaceBody(*body,RigBound.get());
 }

//==============================================================================
/// Configures moving bodies
//==============================================================================
void DSPHChronoLibSC::ConfigMoving(const JChBody *body){
  //Getting moving boundaries
    const double mass=body->GetMass();
    const tdouble3 center=body->GetCenter(); //not used as it is read from a .obj mesh
    ChVector<> COG(center.x,center.y,center.z);
    std::string objname=ChData.GetDataDir()+body->GetModelFile();
    if(!body->GetModelFile().empty()&&!fun::FileExists(objname))fun::Run_ExceptioonFileFun("Model file for moving body is missing.",objname);
    auto MovBound=std::make_shared< ChBodyEasyMesh >(objname.c_str(),1000,false,ChData.GetUseCollision(),CollisionCoef,false,(UseSMC?ChMaterialSurface::SMC:ChMaterialSurface::NSC));
    //if mesh is not read, it still creates a body object
    MphysicalSystem->Add(MovBound);  //Added to chrono scene
    ChMatrix33<> global_csys(1); //-global coord system
    MovBound->SetFrame_COG_to_REF(ChFrame<>(COG,global_csys)); //-seting the COG according to DSPH, not by mesh baricenter
    MovBound->SetNameString(body->IdName);
    MovBound->SetBodyFixed(false);
    MovBound->SetMass(mass);
    //Need to attribute collision parameters				
    //Might be important to define object families to disregard contacts between mechanism parts for example
    if(ChData.GetUseCollision()) ConfigSurfaceBody(*body,MovBound.get());
}

//==============================================================================
/// Loads inertia for bodies.
//==============================================================================
void DSPHChronoLibSC::Config_Inertia(){
  for(unsigned cb=0; cb<ChData.GetBodyCount(); cb++){
    const JChBody* body=ChData.GetBody(cb);
    if(body->Type==JChBody::BD_Floating) { 
      const tmatrix3f inertiaf=ToTMatrix3f(body->GetInertia());
      auto fb=MphysicalSystem->SearchBody(body->IdName.c_str());
      fb->SetInertiaXX(ChVector<>(inertiaf.a11,inertiaf.a22,inertiaf.a33));
      fb->SetInertiaXY(ChVector<>(inertiaf.a12,inertiaf.a23,inertiaf.a13));
    }
  } //all bodies are created in chrono
}

//==============================================================================
/// Compute a given timestep for the full chrono physical system
//==============================================================================
bool DSPHChronoLibSC::RunChrono(double timestep,double dt,bool predictor){
  bool err=false;

  //-Update the system from the beginning until FtPause when it is used
  if(ChData.GetFtPause()>MphysicalSystem->GetChTime())MphysicalSystem->SetChTime(ChData.GetFtPause());

  //-Manages moving objects.
  for(unsigned cb=0;cb<ChData.GetBodyCount()&&!err;cb++)if(ChData.GetBody(cb)->Type==JChBody::BD_Moving){
    JChBodyMoving *body=(JChBodyMoving*)ChData.GetBody(cb);
    auto mb=MphysicalSystem->SearchBody(body->IdName.c_str());
    if(mb!=NULL){
      tdouble3 vlin=TDouble3(0);
      tdouble3 vang=TDouble3(0);
      //-Load input data.
      if(body->GetMotionType()==JChBodyMoving::MV_Simple){
        //-Transforming the translation into linear velocities
        vlin=body->GetMotionSimple()/body->GetMotionDt();
      }
      else if(body->GetMotionType()==JChBodyMoving::MV_Matrix){
        const tmatrix4d mvMatrix=body->GetMotionMatrix();
        //-Now we extract the translation and rotation from the matrix and transform into velocities
        vlin=TDouble3(mvMatrix.a14,mvMatrix.a24,mvMatrix.a34)/body->GetMotionDt(); //-Linear velocity.
        vang=TDouble3(mvMatrix.a11,mvMatrix.a22,mvMatrix.a33)/body->GetMotionDt(); //-Angular velocity.
      }
      else if(body->GetMotionType()!=JChBodyMoving::MV_None)fun::Run_ExceptioonFun(fun::PrintStr("Type of predefined motion applied to object \'%s\' is unknown.",body->IdName.c_str()));
      //-Applying the velocities to the body
      mb->SetPos_dt(ChVector<>(vlin.x,vlin.y,vlin.z));
      mb->SetWvel_loc(ChVector<>(vang.x,vang.y,vang.z));
    }
    else fun::Run_ExceptioonFun(fun::PrintStr("The moving object \'%s\' is missing.",body->IdName.c_str()));
  }


  //-Manages floating objects.
  double yref=0.0;
  for(unsigned cb=0;cb<ChData.GetBodyCount()&&!err; cb++)if(ChData.GetBody(cb)->Type==JChBody::BD_Floating){
    JChBodyFloating *body=(JChBodyFloating*)ChData.GetBody(cb);
    auto fb=MphysicalSystem->SearchBody(body->IdName.c_str());
    //-Do time step and position, velocity e angular velocity of the floating object  
    if(fb!=NULL && body->GetInputData()){
      //-Apply imposed velocity.
      ApplyImposedVel(*body,fb.get());
      //-Apply motion constraints.
      if(!body->GetMotionFree()){
        const tint3 mov=body->GetTranslationFree();
        const ChVector<> vlin=fb->GetPos_dt();
        fb->SetPos_dt(ChVector<>((!mov.x?0:vlin.x()),(!mov.y?0:vlin.y()),(!mov.z?0:vlin.z())));
        const tint3 rot=body->GetRotationFree();
        const ChVector<> vang=fb->GetWvel_loc();
        fb->SetWvel_loc(ChVector<>((!rot.x?0:vang.x()),(!rot.y?0:vang.y()),(!rot.z?0:vang.z())));
      }
      //-Apply forces.
      {
        //-Load input data.
        const tfloat3 face=body->GetInputFace();
        const tfloat3 fomegaace=body->GetInputFomegaAce();
        //-Reset forces
        fb->Empty_forces_accumulators();
        //-Computing force
        ChVector<> force=fb->GetMass()*ChVector<>(face.x,face.y,face.z);
        //-Multiplies the force for the coefficient introduced by user.
        const tfloat3 fcoef=body->GetScaleForce();
        if(fcoef!=TFloat3(FLT_MAX))force=ChVector<>(force.x()*fcoef.x,force.y()*fcoef.y,force.z()*fcoef.z);
        //-Accumulate force
        fb->Accumulate_force(force,fb->GetPos(),false);
        //-Computing torque			
        ChVector<> torque=ChVector<>(fomegaace.x*fb->GetInertiaXX().x()+fb->GetInertiaXY().x()*fomegaace.y+fb->GetInertiaXY().z()*fomegaace.z,
          fomegaace.y*fb->GetInertiaXX().y()+fb->GetInertiaXY().x()*fomegaace.x+fb->GetInertiaXY().y()*fomegaace.z,
          fomegaace.z*fb->GetInertiaXX().z()+fb->GetInertiaXY().z()*fomegaace.x+fb->GetInertiaXY().y()*fomegaace.y);
        fb->Accumulate_torque(torque,false);
        yref=fb->GetPos().y();
      }
    }
    else fun::Run_ExceptioonFun(fun::PrintStr("The floating object \'%s\' is missing.",body->IdName.c_str()));
  }

  //-Establishes the variable coefficients to the link objects.
  if(ChData.GetUseVariableCoeff())SetVariableCoeff();
  
  //-PERFORM SIMULATION UP TO chronoTime  
  MphysicalSystem->DoFrameDynamics(timestep);

  //-Clear motion data of moving bodies (it is important when body stop).
  if(!predictor){
    for(unsigned cb=0;cb<ChData.GetBodyCount()&&!err;cb++)if(ChData.GetBody(cb)->Type==JChBody::BD_Moving)((JChBodyMoving*)ChData.GetBody(cb))->ResetMotion();
  }

  //-Writting data back to the floating bodies in DSPH
  for(unsigned cb=0;cb<ChData.GetBodyCount()&&!err;cb++)if(ChData.GetBody(cb)->Type==JChBody::BD_Floating){
    JChBodyFloating *body=(JChBodyFloating*)ChData.GetBody(cb);
    auto fb=MphysicalSystem->SearchBody(body->IdName.c_str());
    if(Simulate2D){
      ChVector<> Pos=ChVector<>(fb->GetPos().x(),yref,fb->GetPos().z());
      ChVector<> Vel=ChVector<>(fb->GetPos_dt().x(),0.0,fb->GetPos_dt().z());
      ChVector<> Omega=ChVector<>(0.0,fb->GetWvel_par().y(),0.0);
      fb->SetPos(Pos);
      fb->SetPos_dt(Vel);
      fb->SetWvel_par(Omega);
    }
    //-Obtains: center, fvel, fomega.
    const tdouble3 center=TDouble3(fb->GetPos().x(),fb->GetPos().y(),fb->GetPos().z());
    const tfloat3  fvel=ToTFloat3(TDouble3(fb->GetPos_dt().x(),fb->GetPos_dt().y(),fb->GetPos_dt().z()));
    const tfloat3  fomega=ToTFloat3(TDouble3(fb->GetWvel_par().x(),fb->GetWvel_par().y(),fb->GetWvel_par().z()));
    //-Store output data.
    body->SetOutputData(center,fvel,fomega);
  }

  if(!err)RunState=RSTATE_Results;
  return(!err);
}

//==============================================================================
/// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
//==============================================================================
void DSPHChronoLibSC::SaveForcesHead(){
  if(ChData.GetBodyCount()){ //-Creating a file for body reactions to be written on during the run
    ChStreamOutAsciiFile mfileo((DirOut+"ChronoBody_forces.csv").c_str());
    mfileo<<"Time;";
    auto myiter=MphysicalSystem->Get_bodylist().begin();
    while(myiter!=MphysicalSystem->Get_bodylist().end()){
      mfileo<<"Body_"<<(*myiter)->GetName()<<"_fx;fy;fz;mx;my;mz;";
      ++myiter;
    }
    mfileo<<"\n";
  }

  if(ChData.GetLinkCount()){ //-Creating a file for link reactions to be written on during the run
    ChStreamOutAsciiFile mfileo((DirOut+"ChronoLink_forces.csv").c_str());
    mfileo<<"Time;";
    auto myiter=MphysicalSystem->Get_linklist().begin();
    while(myiter!=MphysicalSystem->Get_linklist().end()){
      mfileo<<"Link_"<<(*myiter)->GetName()<<"_fx;fy;fz;mx;my;mz;";
      ++myiter;
    }
    mfileo<<"\n";
  }
}

//==============================================================================
/// Saves forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
//==============================================================================
void DSPHChronoLibSC::SaveForces(){
  //- Writting body forces to file directly
  if(ChData.GetBodyCount()){
    ChStreamOutAsciiFile mfileo((DirOut+"ChronoBody_forces.csv").c_str(),std::ios::app);
    mfileo<<MphysicalSystem->GetChTime()<<";";
    auto myiter=MphysicalSystem->Get_bodylist().begin();
    while(myiter!=MphysicalSystem->Get_bodylist().end()){
      mfileo<<(*myiter)->Get_accumulated_force().x()<<";"<<(*myiter)->Get_accumulated_force().y()<<";"<<(*myiter)->Get_accumulated_force().z()<<";"<<(*myiter)->Get_accumulated_torque().x()<<";"<<(*myiter)->Get_accumulated_torque().y()<<";"<<(*myiter)->Get_accumulated_torque().z()<<";";
      ++myiter;
    }
    mfileo<<"\n";
  }

  //- Writting link forces to file directly
  if(ChData.GetLinkCount()){
    ChStreamOutAsciiFile mfileo((DirOut+"ChronoLink_forces.csv").c_str(),std::ios::app);
    mfileo<<MphysicalSystem->GetChTime()<<";";
    auto myiter=MphysicalSystem->Get_linklist().begin();
    while(myiter!=MphysicalSystem->Get_linklist().end()){
      mfileo<<(*myiter)->Get_react_force().x()<<";"<<(*myiter)->Get_react_force().y()<<";"<<(*myiter)->Get_react_force().z()<<";"<<(*myiter)->Get_react_torque().x()<<";"<<(*myiter)->Get_react_torque().y()<<";"<<(*myiter)->Get_react_torque().z()<<";";
      ++myiter;
    }
    mfileo<<"\n";
  }
}

//==============================================================================
/// Obtains positions of Spring link.
//==============================================================================
bool DSPHChronoLibSC::GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const {
  bool err=true;
  std::shared_ptr<ChLinkBase> link=MphysicalSystem->SearchLink(linkname.c_str());
  if(link){
    //-Linear spring
    if(auto ltype=std::dynamic_pointer_cast<ChLinkSpring>(link)){
      ChVector<> pt1=ltype->GetEndPoint1Abs();
      ChVector<> pt2=ltype->GetEndPoint2Abs();
      p1=TDouble3(pt1.x(),pt1.y(),pt1.z());
      p2=TDouble3(pt2.x(),pt2.y(),pt2.z());
      err=false;
    }
    //-Coulomb damping
    else if(auto ltype=std::dynamic_pointer_cast<ChLinkCoulombDamping>(link)){
      ChVector<> pt1=ltype->GetEndPoint1Abs();
      ChVector<> pt2=ltype->GetEndPoint2Abs();
      p1=TDouble3(pt1.x(),pt1.y(),pt1.z());
      p2=TDouble3(pt2.x(),pt2.y(),pt2.z());
      err=false;
    }
  }
  else fun::Run_ExceptioonFun(fun::PrintStr("The SpringLink \'%s\' is missing.",linkname.c_str()));
  return(err);
}

//==============================================================================
/// Obtains RestLength of Spring link.
//==============================================================================
double DSPHChronoLibSC::GetSpringLinkRestLength(const std::string &linkname)const {
  double ret=0;
  std::shared_ptr<ChLinkBase> link=MphysicalSystem->SearchLink(linkname.c_str());
  if(link){
    if(auto ltype=std::dynamic_pointer_cast<ChLinkSpring>(link)) ret=ltype->Get_SpringRestLength();
    else if(auto ltype=std::dynamic_pointer_cast<ChLinkCoulombDamping>(link)) ret=ltype->Get_SpringRestLength();
  }
  else fun::Run_ExceptioonFun(fun::PrintStr("The SpringLink \'%s\' is missing.",linkname.c_str()));
  return(ret);
}

//==============================================================================
/// Modifies RestLength of Spring link.
//==============================================================================
void DSPHChronoLibSC::SetSpringLinkRestLength(const std::string &linkname,double restlength)const {
  std::shared_ptr<ChLinkBase> link=MphysicalSystem->SearchLink(linkname.c_str());
  if(link){
    if(auto ltype=std::dynamic_pointer_cast<ChLinkSpring>(link)) ltype->Set_SpringRestLength(restlength);
    else if(auto ltype=std::dynamic_pointer_cast<ChLinkCoulombDamping>(link))ltype->Set_SpringRestLength(restlength);
  }
  else fun::Run_ExceptioonFun(fun::PrintStr("The SpringLink \'%s\' is missing.",linkname.c_str()));
}

//==============================================================================
/// Obtains center of body.
//==============================================================================
bool DSPHChronoLibSC::GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const {
  bool err=true;
  auto body=MphysicalSystem->SearchBody(bodyname.c_str());
  if(body){
    pcen=TDouble3(body->GetPos().x(),body->GetPos().y(),body->GetPos().z());
    err=false;
  }
  else fun::Run_ExceptioonFun(fun::PrintStr("The body \'%s\' is missing.",bodyname.c_str()));
  return(err);
}