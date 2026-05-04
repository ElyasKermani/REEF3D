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

#include"contact_force.h"
#include<cmath>

contact_force::~contact_force()
{
}

double contact_force::effective_young_modulus(double E1, double E2, double nu1, double nu2)
{
    return 1.0 / ((1.0 - nu1*nu1)/E1 + (1.0 - nu2*nu2)/E2);
}

double contact_force::effective_radius(double R1, double R2)
{
    return 1.0 / (1.0/R1 + 1.0/R2);
}

double contact_force::hertz_stiffness(double E_eff, double R_eff)
{
    return E_eff * std::sqrt(R_eff);
}
