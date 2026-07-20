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

Non-linear Hertzian elastic contact with LIGGGHTS-style normal damping
and Mindlin tangential stiffness/damping.
--------------------------------------------------------------------*/

#ifndef CONTACT_FORCE_HERTZ_H_
#define CONTACT_FORCE_HERTZ_H_

#include"contact_force.h"

class lexer;

using namespace std;

class contact_force_hertz : public contact_force
{
public:

    contact_force_hertz(lexer*);
    virtual ~contact_force_hertz();

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

    double E;     // Young modulus (Pa); normal contact stiffness follows Hertz
    double nu;    // Poisson ratio (nondim.); used with E for effective contact moduli
    double mu;    // Coulomb friction coefficient (nondim.)
    double cor;   // normal coefficient of restitution (nondim.)
};

#endif
