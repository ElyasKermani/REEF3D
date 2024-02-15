/*
 <DUALSPHYSICS>  Copyright (c) 2019, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics.

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

/// \file DSPHChronoLib.h \brief Declares the class \ref DSPHChronoLib which is the interface between DualSPHysics and Chrono.

#ifndef DSPHCHRONOLIB_H
#define DSPHCHRONOLIB_H

#include "TypesDef.h"
#include "JChronoData.h"
#include <string>
#include <iostream>
#include <memory>

//-Forward declarations to avoid including chrono classes.
namespace chrono {
  class ChSystem;
  class ChMaterialSurface;
  class ChBody;
};

//##############################################################################
//# DSPHChronoLib
//##############################################################################
/// \brief Defines the class interface between DualSPHysics and Chrono.
class DSPHChronoLib {
public:
  //-States of execution.
  typedef enum { RSTATE_Init,RSTATE_Loading,RSTATE_Results }TpRunState;
  const std::string version;       ///<DualSPHysics version
  const std::string DsphChVersion; ///<Interface version

protected:
  const std::string ClassName;

  //-Chrono physical system.
  std::string DirOut;
  bool Simulate2D;          ///<True for 2D Simulations.
  bool UseOmp;              ///<Indicates if use of ChronoEngine_Multicore module is enabled.
  int OmpThreads;           ///<Threads number used by OpenMP.
  unsigned SolverIx;        ///<Indicates the index of chrono solver enum.
  unsigned TimeStepperIx;   ///<Indicates the index of chrono timestepper enum.
  bool UseSMC;              ///<True if it is using SMC (SMooth Contacts) 
  JChronoData ChData;
  TpRunState RunState;
  double CollisionCoef;

  /// Constructor
  DSPHChronoLib(const JChronoData &chdata);

  /// Initialisation of variables
  void Reset();

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  virtual void SaveForcesHead(){};
  
  /// Establishes the variable coefficients to the link objects.
  virtual void SetVariableCoeff(){};
  
  /// Adds the material properties to a object to enable collisions
  void ConfigSurfaceBody(const JChBody &body,chrono::ChBody *chbody);

  /// Adds the initial velocity.
  void ApplyInitialVel(const JChBody &body,chrono::ChBody *chbody);

  /// Adds the imposed velocity.
  void ApplyImposedVel(const JChBodyFloating &body,chrono::ChBody *chbody);

public:

  /// Loads data for bodies and configures objects.
  virtual void Config(std::string dirout,bool svdata,bool simulate2d){};

  /// Loads inertia for bodies.
  virtual void Config_Inertia(){};

  /// Compute a single timestep for each floating and moving body.
  virtual bool RunChrono(double timestep,double dt,bool predictor)=0;
 
  /// Saves forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  virtual void SaveForces(){};

  /// Obtains positions of Spring link.
  virtual bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const=0;

  /// Obtains RestLength of Spring link.
  virtual double GetSpringLinkRestLength(const std::string &linkname)const=0;

  /// Modifies RestLength of Spring link.
  virtual void SetSpringLinkRestLength(const std::string &linkname,double restlength)const{};

  /// Obtains center of body.
  virtual bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const=0;

  /// Returns pointer to ChronoData object.
  const JChronoData* GetChronoData(){ return(&ChData); }

  /// Loads floating data to calculate coupling with Chrono.
  bool SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace);

  /// Loads imposed velocity for floating to calculate coupling with Chrono.
  bool SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang); 
  
  /// Obtains floating data from coupling with Chrono.
  bool GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const;

  /// Loads motion data to calculate coupling with Chrono.
  bool SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt);
};

//##############################################################################
//# DSPHChronoLibSC
//##############################################################################
/// \brief Defines the class for single-core executions.
class DSPHChronoLibSC : public DSPHChronoLib {
private:
  chrono::ChSystem *MphysicalSystem; ///<Pointer to Chrono System

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  void SaveForcesHead();
  
  /// Configures floating bodies
  void ConfigFloating(const JChBody* body);

  /// Configures moving bodies
  void ConfigMoving(const JChBody* body);

  /// Configures fixed bodies
  void ConfigFixed(const JChBody* body);

public:
  /// Constructor
  DSPHChronoLibSC(const JChronoData &chdata);

  /// Destructor
  ~DSPHChronoLibSC();

  /// Loads data for bodies and configures objects.
  void Config(std::string dirout,bool svdata,bool simulate2d);

  /// Loads inertia for bodies.
  void Config_Inertia();

  /// Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);

  /// Saves forces for each body and link (ChronoLink_forces.csv,ChronoBody_forces.csv).
  void SaveForces();

  /// Obtains positions of Spring link.
  bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const;

  /// Obtains RestLength of Spring link.
  double GetSpringLinkRestLength(const std::string &linkname)const;

  /// Modifies RestLength of Spring link.
  void SetSpringLinkRestLength(const std::string &linkname,double restlength)const;

  /// Obtains center of body.
  bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const;
};
#endif //!DSPHCHRONOLIB_H