/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef CONTACT_FORCE_LINEAR_H_
#define CONTACT_FORCE_LINEAR_H_

#include"contact_force.h"

class lexer;

using namespace std;

class contact_force_linear : public contact_force
{
public:

    contact_force_linear(lexer*);
    virtual ~contact_force_linear();

    void compute(lexer*, ghostcell*,
                 sixdof_obj *obj1, sixdof_obj *obj2,
                 const Eigen::Vector3d &contact_point,
                 const Eigen::Vector3d &normal,
                 double overlap,
                 ContactHistory &history,
                 double dt_contact,
                 bool finalize,
                 Eigen::Vector3d &force, Eigen::Vector3d &torque) override;

private:

    double kn;       // normal spring stiffness (N/m)
    double kt;       // tangential spring stiffness (N/m)
    double cn;       // normal viscous damping (N·s/m); overridden if use_cor
    double ct;       // tangential viscous damping (N·s/m); overridden if use_cor
    double mu;       // Coulomb friction coefficient (nondim.)
    double en;       // normal coefficient of restitution (nondim.)
    double et;       // tangential restitution (nondim.); -1 means use en for tangential damping
    bool use_cor;    // if true, derive cn/ct from en/et and effective mass
};

#endif
