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

#ifndef _6DOF_BOUNDARY_WALL_H_
#define _6DOF_BOUNDARY_WALL_H_

#include<Eigen/Dense>

class lexer;
class ghostcell;

class boundary_wall_6DOF
{
public:
    boundary_wall_6DOF(int wall_id, const Eigen::Vector3d &normal, const Eigen::Vector3d &point_on_wall);
    
    // Wall identifier (e.g., 0=xmin, 1=xmax, 2=ymin, etc.)
    int id;
    
    // Wall normal vector (pointing inward to the domain)
    Eigen::Vector3d normal_vector;
    
    // A point on the wall
    Eigen::Vector3d reference_point;
    
    // Physical properties
    double friction_coefficient;
    double restitution_coefficient;
    double stiffness;
    double damping;
};

#endif 