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

#ifndef CHLINKCOULOMBDAMPING_H
#define CHLINKCOULOMBDAMPING_H

#include "chrono/physics/ChLinkSpring.h"

namespace chrono {

/// Class for spring-damper systems, acting along the polar
/// distance of two markers

class ChApi ChLinkCoulombDamping : public ChLinkSpring {
  protected:
    double CoulombDamping; ///< Coulomb damping value (it is not applied when it is zero).  //<vs_dsphchanges>

  public:
    ChLinkCoulombDamping();
    ChLinkCoulombDamping(const ChLinkCoulombDamping& other);
    virtual ~ChLinkCoulombDamping(){}

    /// "Virtual" copy constructor (covariant return type).
    virtual ChLinkCoulombDamping* Clone() const override { return new ChLinkCoulombDamping(*this); }

    /// Specialized initialization for springs, given the two bodies to be connected,
    /// the positions of the two anchor endpoints of the spring (each expressed
    /// in body or abs. coordinates) and the imposed rest length of the spring.
    /// NOTE! As in ChLinkMarkers::Initialize(), the two markers are automatically
    /// created and placed inside the two connected bodies.
    void Initialize(
        std::shared_ptr<ChBody> mbody1,  ///< first body to link
        std::shared_ptr<ChBody> mbody2,  ///< second body to link
        bool pos_are_relative,           ///< true: following pos. are relative to bodies
        ChVector<> mpos1,                ///< position of spring endpoint, for 1st body (rel. or abs., see flag above)
        ChVector<> mpos2,                ///< position of spring endpoint, for 2nd body (rel. or abs., see flag above)
        bool auto_rest_length = true,    ///< if true, initializes the rest-length as the distance between mpos1 and mpos2
        double mrest_length = 0,         ///< imposed rest_length (no need to define, if auto_rest_length=true.)
        double coulombdamping = 0        ///< Coulomb damping value (it is not applied when it is zero). <vs_dsphchanges>
        );

    void UpdateForces(double mytime);

  };


CH_CLASS_VERSION(ChLinkCoulombDamping,0)

}  // end namespace chrono

#endif
