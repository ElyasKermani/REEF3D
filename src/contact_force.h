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

#ifndef CONTACT_FORCE_H_
#define CONTACT_FORCE_H_

#include<Eigen/Dense>

class lexer;
class ghostcell;
class sixdof_obj;
struct ContactHistory;

using namespace std;

class contact_force
{
public:

    virtual ~contact_force();

    virtual void compute(lexer*, ghostcell*,
                         sixdof_obj *obj1, sixdof_obj *obj2,
                         const Eigen::Vector3d &contact_point,
                         const Eigen::Vector3d &normal,
                         double overlap,
                         ContactHistory &history,
                         Eigen::Vector3d &force, Eigen::Vector3d &torque)=0;

protected:

    // Ei: Young modulus (Pa); nui: Poisson ratio (nondim.) for bodies 1 and 2
    static double effective_young_modulus(double E1, double E2, double nu1, double nu2);
    // Equivalent sphere radii R1,R2 at contact (m)
    static double effective_radius(double R1, double R2);
    // Curved-contact stiffness scale from composite E_eff (Pa) and contact radius R_eff (m)
    static double hertz_stiffness(double E_eff, double R_eff);
};

#endif
