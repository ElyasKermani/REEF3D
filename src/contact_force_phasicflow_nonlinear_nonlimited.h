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

phasicFlow Hertz-Mindlin model with no tangential damping. Tangential
history is NOT rescaled when the Coulomb friction limit is reached.
--------------------------------------------------------------------*/

#ifndef CONTACT_FORCE_PHASICFLOW_NONLINEAR_NONLIMITED_H_
#define CONTACT_FORCE_PHASICFLOW_NONLINEAR_NONLIMITED_H_

#include"contact_force.h"

class lexer;

using namespace std;

class contact_force_phasicflow_nonlinear_nonlimited : public contact_force
{
public:

    contact_force_phasicflow_nonlinear_nonlimited(lexer*);
    virtual ~contact_force_phasicflow_nonlinear_nonlimited();

    void compute(lexer*, ghostcell*,
                 sixdof_obj *obj1, sixdof_obj *obj2,
                 const Eigen::Vector3d &contact_point,
                 const Eigen::Vector3d &normal,
                 double overlap,
                 ContactHistory &history,
                 Eigen::Vector3d &force, Eigen::Vector3d &torque) override;

private:

    double E;     // Young modulus (Pa); Hertz normal stiffness
    double nu;    // Poisson ratio (nondim.); shear modulus G = E/(1+nu) in this kernel
    double mu;    // Coulomb friction coefficient (nondim.)
    double cor;   // normal coefficient of restitution (nondim.)
};

#endif
