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

#ifndef WALL_CONTACT_H_
#define WALL_CONTACT_H_

#include"boundary_contact.h"

class wall_contact : public boundary_contact
{
public:

    wall_contact();
    virtual ~wall_contact();

    void apply(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj) override;

private:

    double Kw_;   // wall stiffness scale [N/m^2]; effective k_eff = Kw_ * radius
    double mu_;   // Coulomb friction coefficient (nondim.)
    double xi_;   // damping ratio (nondim.)
    double nu_;   // friction velocity smoothing scale [m/s]
};

#endif
