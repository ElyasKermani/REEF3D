/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef HEAT_MATERIAL_H_
#define HEAT_MATERIAL_H_

// Temperature dependent properties of water and air (temperature in deg C), shared by the
// heat-coupled property updates (H 3). Piecewise linear interpolation of tabulated values,
// constant beyond the ends of the tables.

class lexer;

namespace heat_material
{
    double water_density(double);       // kg/m3
    double air_density(double);         // kg/m3
    double water_viscosity(double);     // kinematic, m2/s
    double air_viscosity(double);       // kinematic, m2/s

    double water_beta(double);          // thermal expansion coefficient -(1/rho) drho/dT, 1/K

    // density and kinematic viscosity of level-set phase 1 (phi > 0) or 2 at temperature T, following
    // H 3 (1: tables, 2: linear equation of state rho0 (1 - beta (T - T0)), constant viscosity) and H 9
    // (which phase is water). H 3 = 2: rho0 = W 1 / W 3, T0 = H 50 of the phase, beta from H 4 or by
    // default the tabulated value for water at T0 and 1/T0 (ideal gas, T0 in K) for air.
    double density(lexer*, int, double);
    double viscosity(lexer*, int, double);
}

#endif
