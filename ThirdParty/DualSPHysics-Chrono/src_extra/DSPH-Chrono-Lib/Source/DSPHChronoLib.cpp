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

/// \file DSPHChronoLib.cpp \brief Implements the class \ref DSPHChronoLib.

#include "DSPHChronoLib.h"
#include "Functions.h"
#include "FunDSPHChrono.h"

#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

//==============================================================================
//-------- Physics dependencies ---------
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChMaterialSurfaceNSC.h"
#include "chrono/solver/ChSolver.h"


//-------- Irrlicht dependencies ---------
#ifdef _IRRLICHT_MODULE //Irrlicht module to visualize the chrono objects
#include "chrono_irrlicht/ChIrrApp.h" 
#endif
//-------- Core dependencies ---------
#include "chrono/core/ChVector.h"

#include <cmath>
#include <string>
#include <fstream>
#include <memory.h>

//From chrono demos
#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

//-Use the namespace of Chrono.
using namespace ::chrono;
using namespace ::chrono::collision;

namespace chrono {
  ChApi const ChVector<double> VNULL(0.,0.,0.);
  ChApi const ChVector<double> VECT_X(1.,0.,0.);
  ChApi const ChVector<double> VECT_Y(0.,1.,0.);
  ChApi const ChVector<double> VECT_Z(0.,0.,1.);
  ChApi const ChQuaternion<double> QNULL(0.,0.,0.,0.);
  ChApi const ChQuaternion<double> QUNIT(1.,0.,0.,0.);
  ChApi const ChCoordsys<double> CSYSNULL(VNULL, QNULL);
  ChApi const ChCoordsys<double> CSYSNORM(VNULL, QUNIT);
}

//-------------------------------------------------------------------------------------
/// ------------------------------DSPHChronoLib----------------------------------------
//-------------------------------------------------------------------------------------

//==============================================================================
/// Constructor
//==============================================================================
DSPHChronoLib::DSPHChronoLib(const JChronoData &chdata)
  : version("v5.0.231"), ClassName("DSPHChronoLib"),DsphChVersion("4.21") {
	Reset();
  //-Set path to Chrono data directory.
  SetChronoDataPath(CHRONO_DATA_DIR);
  ChData=chdata;
  OmpThreads=chdata.GetOmpThreads();
  UseSMC=chdata.GetContactMethod()==JChronoData::SMC?true:false;
  UseOmp=(OmpThreads<=1?false:true);
  //-Config JChronoData
  ChData.SetExecMode(UseOmp?"Multi Core":"Single Core");
  //-Set index solver
  SolverIx=fdsch::GetSolverChIxAuto(UseSMC,ChData);
}

//==============================================================================
/// Initialisation of variables
//==============================================================================
void DSPHChronoLib::Reset(){
  UseSMC=UseOmp=false;
  OmpThreads=1;
  CollisionCoef=0;
  SolverIx=0;
}

//==============================================================================
/// Loads floating data to calculate coupling with Chrono.
//==============================================================================
bool DSPHChronoLib::SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace){
  if(RunState!=RSTATE_Loading){
    for(unsigned cb=0; cb<ChData.GetBodyCount(); cb++)if(ChData.GetBody(cb)->Type==JChBody::BD_Floating){
      ((JChBodyFloating *)ChData.GetBody(cb))->ResetInputData();
    }
    RunState=RSTATE_Loading;
  }
  JChBodyFloating *body=(JChBodyFloating*)ChData.GetBodyFloating(mkbound);
  if(body)body->SetInputData(face,fomegaace);
  return(body!=NULL);
}

//==============================================================================
/// Loads imposed velocity for floating to calculate coupling with Chrono.
//==============================================================================
bool DSPHChronoLib::SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang){
  JChBodyFloating *body=(JChBodyFloating*)ChData.GetBodyFloating(mkbound);
  if(body)body->SetInputDataVel(vlin,vang);
  return(body!=NULL);
}

//==============================================================================
/// Obtains floating data from coupling with Chrono.
//==============================================================================
bool DSPHChronoLib::GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const {
  const JChBodyFloating *body=(RunState==RSTATE_Results?ChData.GetBodyFloating(mkbound):NULL);
  if(body){
    fcenter=body->GetOutputCenter();
    fvel=body->GetOutputVel();
    fomega=body->GetOutputOmega();
  }
  return(body!=NULL);
}

//==============================================================================
/// Loads motion data to calculate coupling with Chrono.
//==============================================================================
bool DSPHChronoLib::SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt){
  JChBodyMoving *body=(JChBodyMoving*)ChData.GetBodyMoving(mkbound);
  if(body){
    if(simple)body->SetMotionSimple(stepdt,msimple);
    else body->SetMotionMatrix(stepdt,mmatrix);
  }
  return(body!=NULL);
}

//==============================================================================
  /// Adds the material properties to a object to enable collisions
//==============================================================================
void DSPHChronoLib::ConfigSurfaceBody(const JChBody &body,ChBody *chbody){
  //Need to define material properties
  if(UseSMC){ //SMooth Contacts
    auto mmaterial=std::make_shared<ChMaterialSurfaceSMC>();
    mmaterial->SetYoungModulus(body.GetYoung());
    mmaterial->SetPoissonRatio(body.GetPoisson());
    mmaterial->SetFriction(body.GetKfric());
    mmaterial->SetRestitution(body.GetRestitu());
    mmaterial->SetImposeFric(body.GetImposeFric());
    mmaterial->SetAdhesion(0);
    chbody->SetMaterialSurface(mmaterial);
  }
  else {//Non Smooth Contacts
    auto mmaterial=std::make_shared<ChMaterialSurfaceNSC>();
    mmaterial->SetFriction(body.GetKfric());
    mmaterial->SetRestitution(body.GetRestitu());
    mmaterial->SetImposeFric(body.GetImposeFric());
    chbody->SetMaterialSurface(mmaterial);
  }
}

//==============================================================================
/// Adds the initial velocity.
//==============================================================================
void DSPHChronoLib::ApplyInitialVel(const JChBody &body,ChBody *chbody){
  //-Loads initial velocities
  const tfloat3 vlinini=body.GetLinearVelini();
  const tfloat3 vangini=body.GetAngularVelini();
  //-Apply initial linear velocity
  if(vlinini!=TFloat3(0)) chbody->SetPos_dt(ChVector<>(vlinini.x,vlinini.y,vlinini.z));
  //-Apply initial angular velocity
  if(vangini!=TFloat3(0)) chbody->SetWvel_loc(ChVector<>(vangini.x,vangini.y,vangini.z));
}

//==============================================================================
/// Adds the imposed velocity.
//==============================================================================
void DSPHChronoLib::ApplyImposedVel(const JChBodyFloating &body,ChBody *chbody){
  //-Loads velocities
  const tfloat3 vlin=body.GetInputLinearVel();
  const tfloat3 vang=body.GetInputAngularVel();
  //-Apply linear velocity
  if(vlin!=TFloat3(FLT_MAX)){
    const ChVector<> vlin0=chbody->GetPos_dt();
    chbody->SetPos_dt(ChVector<>((vlin.x!=FLT_MAX?vlin.x:vlin0.x()),(vlin.y!=FLT_MAX?vlin.y:vlin0.y()),(vlin.z!=FLT_MAX?vlin.z:vlin0.z())));
  }
  //-Apply angular velocity
  if(vang!=TFloat3(FLT_MAX)){
    const ChVector<> vang0=chbody->GetWvel_loc();
    chbody->SetWvel_loc(ChVector<>((vang.x!=FLT_MAX?vang.x:vang0.x()),(vang.y!=FLT_MAX?vang.y:vang0.y()),(vang.z!=FLT_MAX?vang.z:vang0.z())));
  }
}