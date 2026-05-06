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

#ifndef GROUND_CONTACT_H_
#define GROUND_CONTACT_H_

#include"boundary_contact.h"

class ground_contact : public boundary_contact
{
public:

    ground_contact();
    virtual ~ground_contact();

    void apply(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj) override;

private:

    double zg_;   // seabed z (inertial); penalty plane
    double Kg_;   // stiffness scale [N/m^2]; effective k_eff = Kg_ * radius
    double mu_;   // Coulomb friction coefficient (nondim.)
    double xi_;   // damping ratio for normal spring (nondim.)
    double nu_;   // friction velocity smoothing scale [m/s]
};

#endif
