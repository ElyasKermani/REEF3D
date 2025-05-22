/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include "6DOF_boundary_wall.h"
#include "lexer.h"
#include "ghostcell.h"

boundary_wall_6DOF::boundary_wall_6DOF(int wall_id, const Eigen::Vector3d &normal, const Eigen::Vector3d &point_on_wall)
: id(wall_id), normal_vector(normal), reference_point(point_on_wall)
{
    // Default physical properties (can be set from configuration)
    friction_coefficient = 0.3;
    restitution_coefficient = 0.7;
    stiffness = 1.0e6;
    damping = 1.0e4;
} 