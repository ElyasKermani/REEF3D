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

/// \file FunDSPHChrono.cpp \brief Implements geometry functions for DSPH-Chrono.

#include "FunDSPHChrono.h"
#include "Functions.h"

#include <cstdio>
#include <climits>
#include <algorithm>

namespace fdsch{

  //==============================================================================
  // Convierte un cuaternion a angulos de euler.
  // Converts a quaternion in euler-angles.
  //==========================================================================
  tdouble3 QToEuler(chrono::ChQuaternion<> q){
    tdouble3 euler;
    const double sq0=q.e0()*q.e0();
    const double sq1=q.e1()*q.e1();
    const double sq2=q.e2()*q.e2();
    const double sq3=q.e3()*q.e3();

    euler.x=+atan2(2*(q.e2()*q.e3()+q.e0()*q.e1()),sq3-sq2-sq1+sq0);//- roll
    euler.y=-asin (2*(q.e1()*q.e3()-q.e0()*q.e2()));                //- pitch
    euler.z=+atan2(2*(q.e1()*q.e2()+q.e3()*q.e0()),sq1+sq0-sq3-sq2);//- yaw
  
    if(std::isnan(euler.y))euler.y=(-90*PI/180);//-Avoid the limitation of euler-angles

    return (euler);
  }

  //==============================================================================
  // Comprueba si hay Not-A-Number en un array de tdouble3.
  // Check if there are Not-A-Number in a tdouble3 array.
  //==============================================================================
  bool CheckForNaNs(const tdouble3 *v,const unsigned size){
    bool err=false;
    for(unsigned c=0;c<size && !err;c++){
      tdouble3 val=v[c];
      if(std::isnan(val.x) || std::isnan(val.y) || std::isnan(val.z))
        err=true;
    }
    return (err);
  }
  
  //==============================================================================
  // Comprueba si hay Not-A-Number en un array de tfloat3.
  // Check if there are Not-A-Number in a tfloat3 array.
  //==============================================================================
  bool CheckForNaNs(const tfloat3 *v,const unsigned size){
    bool err=false;
    for(unsigned c=0;c<size && !err;c++){
      tfloat3 val=v[c];
      if(std::isnan(val.x) || std::isnan(val.y) || std::isnan(val.z))
        err=true;
    }
    return (err);
  }

  //==============================================================================
  /// Return an array of tdouble3 with the requiered size.
  //==============================================================================
  tdouble3 * GetDouble3Array(const unsigned size){
    tdouble3 *arr=new tdouble3[size];
    memset(arr,0,sizeof(tdouble3)*size);
    return arr;
  }
}



