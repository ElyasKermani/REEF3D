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
// Authors: Ivan Martinez-Estevez
// =============================================================================

#include "chrono/physics/ChLinkCoulombDamping.h"

namespace chrono {

// Register into the object factory, to enable run-time
// dynamic creation and persistence
CH_FACTORY_REGISTER(ChLinkCoulombDamping)

ChLinkCoulombDamping::ChLinkCoulombDamping():ChLinkSpring() {
    CoulombDamping = 0.0; //<vs_dsphchanges>
}

ChLinkCoulombDamping::ChLinkCoulombDamping(const ChLinkCoulombDamping& other) : ChLinkSpring(other) {
    CoulombDamping = other.CoulombDamping; //<vs_dsphchanges>
}

void ChLinkCoulombDamping::Initialize(std::shared_ptr<ChBody> mbody1,
                              std::shared_ptr<ChBody> mbody2,
                              bool pos_are_relative,
                              ChVector<> mpos1,
                              ChVector<> mpos2,
                              bool auto_rest_length,
                              double mrest_length,
                              double coulombdamping) {
    // First, initialize as all constraint with markers.
    // In this case, create the two markers also!.
    ChLinkSpring::Initialize(mbody1,mbody2,pos_are_relative,mpos1,mpos2,auto_rest_length,mrest_length);
    CoulombDamping=coulombdamping; //<vs_dsphchanges>
}


void ChLinkCoulombDamping::UpdateForces(double mytime) {
  // Inherit force computation:
  // also base class can add its own forces.
  ChLinkMarkers::UpdateForces(mytime);

  spr_react = 0.0;
  Vector m_force;
  double deform = Get_SpringDeform();

//<vs_dsphchanges_ini> 
  //Original Chrono code.
  //spr_react = spr_f * mod_f_time->Get_y(ChTime);
  //spr_react -= (spr_k * mod_k_d->Get_y(deform) * mod_k_speed->Get_y(dist_dt)) * (deform);
  //spr_react -= (spr_r * mod_r_d->Get_y(deform) * mod_r_speed->Get_y(dist_dt)) * (dist_dt);

  //Modified Chrono code.
  spr_react = spr_f * mod_f_time->Get_y(ChTime);
    spr_react -= (mod_r_d->Get_y(deform) * mod_r_speed->Get_y(dist_dt)) * (dist_dt>=0? CoulombDamping: -CoulombDamping); // This is enough for link_coulombdamping

  m_force = Vmul(Vnorm(relM.pos), spr_react);

  C_force = Vadd(C_force, m_force);
}

}  // end namespace chrono
