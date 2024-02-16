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

//:#############################################################################
//:# Descripcion:
//:# =============
//:# Conjunto de funciones tipicas para utilizar en el acpolamiento DSPH-Chrono.
//:#
//:# Cambios:
//:# =========
//:# - Implementacion del namespace. (29-06-2021)
//:#############################################################################

/// \file FunDSPHChrono.h \brief Declares functions for DSPH-Chrono.

#ifndef _FunDSPHChrono_
#define _FunDSPHChrono_

#include "TypesDef.h"
#include "JChronoData.h"

#include "chrono/core/ChVector2.h"
#include "chrono/core/ChQuaternion.h"
#include "chrono/solver/ChSolver.h"
#include "chrono/timestepper/ChTimestepper.h"

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <memory>


/// Implements a set of functions for DSPH-Chrono.
namespace fdsch{
  tdouble3 QToEuler(chrono::ChQuaternion<> q);
  bool CheckForNaNs(const tdouble3 *v,const unsigned size);
  bool CheckForNaNs(const tfloat3 *v,const unsigned size);
  tdouble3 *GetDouble3Array(const unsigned size);
  //==============================================================================
  /// Devuelve tdouble3 a partir de un ChVector<>.
  /// Returns a tdouble3 from a ChVector<>.
  //==============================================================================
  inline tdouble3 ChVecToDouble3(chrono::ChVector<> v){
    return TDouble3(v.x(),v.y(),v.z());
  }
  
  //==============================================================================
  /// Devuelve tfloat3 a partir de un ChVector<>.
  /// Returns a tfloat3 from a ChVector<>.
  //==============================================================================
  inline tfloat3 ChVecToFloat3(chrono::ChVector<> v){
    return TFloat3((float)v.x(),(float)v.y(),(float)v.z());
  }

  //==============================================================================
  /// Devuelve un ChVector<> a partir de un tdouble3.
  /// Returns a ChVector<> from a tdouble3.
  //==============================================================================
  inline chrono::ChVector<> Real3ToChVec(tdouble3 v){
    return chrono::ChVector<>(v.x,v.y,v.z);
  }

  //==============================================================================
  /// Devuelve un ChVector<> a partir de un tfloat3.
  /// Returns a ChVector<> from a tfloat3.
  //=============================================================================
  inline chrono::ChVector<> Real3ToChVec(tfloat3 v){
    return chrono::ChVector<>(v.x,v.y,v.z);
  }

  //==============================================================================
  /// Devuelve vector unitario valido del vector or (0,0,0).
  /// Returns a valid unit vector of the vector or (0,0,0).
  //==============================================================================
  inline tdouble3 VecUnitary(const tdouble3 &v){
    const double m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
    return(m? TDouble3(v.x/m,v.y/m,v.z/m): v);
  }

  //==============================================================================
  /// Explicit conversion of scoped enumeration to int (e.g. for streaming).
  //==============================================================================
  #ifdef DISABLE_CHRONO_OMP
  template <typename Enumeration>
  inline auto as_integer(Enumeration const value) -> typename std::underlying_type<Enumeration>::type {
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
  }
  #endif

  ////==============================================================================
  /////Returns the solver index in function of the execution mode
  ////==============================================================================
  inline int GetSolverChIxAuto(bool smc,JChronoData &chdata){
    int index=-1;
    if(smc){
      index=fdsch::as_integer(chrono::ChSolver::Type::MINRES);
      chdata.SetSolver(JChronoData::TpSolverType::MINRES);
    }
    else{
      index=fdsch::as_integer(chrono::ChSolver::Type::BARZILAIBORWEIN);
      chdata.SetSolver(JChronoData::TpSolverType::BB);
    }
    return (index);
  }
}
#endif


