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

#ifndef SIXDOF_COLLISION_H_
#define SIXDOF_COLLISION_H_

#include<Eigen/Dense>
#include<vector>

class lexer;
class ghostcell;
class sixdof_obj;

using namespace std;

class sixdof_collision
{
public:

    sixdof_collision(lexer *p, ghostcell *pgc);
    virtual ~sixdof_collision();
    
    // Calculate collision forces between all 6DOF objects
    void calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
private:

    // Detect collision between two 6DOF objects
    bool detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2, 
                         Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);
    
    // Calculate linear contact force using the linear spring-dashpot model
    void calculate_linear_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                      const Eigen::Vector3d &contact_point, 
                                      const Eigen::Vector3d &normal, 
                                      const double overlap,
                                      Eigen::Vector3d &force, 
                                      Eigen::Vector3d &torque);
    
    // Parameters for the linear collision model
    double spring_constant;    // Normal spring stiffness [N/m]
    double damping_constant;   // Normal damping coefficient [N·s/m]
    double friction_coefficient; // Tangential friction coefficient
    double restitution_coefficient; // Coefficient of restitution
    
    // Simplified bounding sphere collision detection
    double calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2);
    
    // For simplified spherical collision (using bounding spheres)
    vector<double> bounding_radius;
};

#endif 