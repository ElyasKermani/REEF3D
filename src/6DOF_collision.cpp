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
Author: 
--------------------------------------------------------------------*/

#include"6DOF_collision.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>
#include<iostream>

sixdof_collision::sixdof_collision(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF Collision Model startup..."<<endl;
    
    // Initialize collision model parameters from control file parameters
    // These should be added to the parameter file in a real implementation
    // For now, we set default values
    spring_constant = 1.0e6;            // Default stiffness [N/m]
    damping_constant = 1.0e4;           // Default damping [N·s/m]
    friction_coefficient = 0.3;         // Default friction coefficient
    restitution_coefficient = 0.7;      // Default restitution coefficient
    
    // Allocate space for bounding radius for each 6DOF object
    bounding_radius.resize(p->X20);
    
    // Set initial bounding radius for each object
    // In a real implementation, this would be calculated based on the object geometry
    for(int i=0; i<p->X20; ++i)
    {
        bounding_radius[i] = 1.0; // Default radius (should be calculated from geometry)
    }
}

sixdof_collision::~sixdof_collision()
{
}

void sixdof_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Check for collisions between all pairs of objects
    for(int i=0; i<p->X20-1; ++i)
    {
        for(int j=i+1; j<p->X20; ++j)
        {
            // Variables to store collision information
            Eigen::Vector3d contact_point, normal;
            double overlap = 0.0;
            
            // Detect if collision occurred
            if(detect_collision(p, pgc, fb_obj[i], fb_obj[j], contact_point, normal, overlap))
            {
                // Calculate forces and torques from collision
                Eigen::Vector3d force, torque;
                
                // Apply linear contact model
                calculate_linear_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                              contact_point, normal, overlap, 
                                              force, torque);
                
                // Add collision forces to object i (action)
                fb_obj[i]->Xext += force(0);
                fb_obj[i]->Yext += force(1);
                fb_obj[i]->Zext += force(2);
                
                fb_obj[i]->Kext += torque(0);
                fb_obj[i]->Mext += torque(1);
                fb_obj[i]->Next += torque(2);
                
                // Add collision forces to object j (reaction)
                fb_obj[j]->Xext -= force(0);
                fb_obj[j]->Yext -= force(1);
                fb_obj[j]->Zext -= force(2);
                
                fb_obj[j]->Kext -= torque(0);
                fb_obj[j]->Mext -= torque(1);
                fb_obj[j]->Next -= torque(2);
                
                if(p->mpirank==0 && p->count%p->P12==0)
                {
                    cout<<"6DOF Collision detected between objects "<<i<<" and "<<j<<endl;
                    cout<<"  Overlap: "<<overlap<<endl;
                    cout<<"  Force: ["<<force(0)<<", "<<force(1)<<", "<<force(2)<<"]"<<endl;
                }
            }
        }
    }
}

bool sixdof_collision::detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                      Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Using simplified spherical collision detection for this example
    // In a real implementation, this would need to be replaced with proper geometry-based collision detection
    
    int id1 = obj1->n6DOF;
    int id2 = obj2->n6DOF;
    
    // Get object centers
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Calculate distance between centers
    Eigen::Vector3d center_diff = center2 - center1;
    double distance = center_diff.norm();
    
    // Calculate sum of bounding radii
    double sum_radii = bounding_radius[id1] + bounding_radius[id2];
    
    // Check if objects are overlapping
    if(distance < sum_radii)
    {
        // Calculate overlap
        overlap = sum_radii - distance;
        
        // Calculate contact normal
        if(distance > 1.0e-10)
        {
            normal = center_diff / distance;
        }
        else
        {
            // If centers are too close, use a default normal
            normal = Eigen::Vector3d(0.0, 0.0, 1.0);
        }
        
        // Calculate contact point (midpoint of overlap)
        contact_point = center1 + normal * (bounding_radius[id1] - 0.5 * overlap);
        
        return true;
    }
    
    return false;
}

void sixdof_collision::calculate_linear_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                  const Eigen::Vector3d &contact_point, 
                                                  const Eigen::Vector3d &normal, 
                                                  const double overlap,
                                                  Eigen::Vector3d &force, 
                                                  Eigen::Vector3d &torque)
{
    int id1 = obj1->n6DOF;
    int id2 = obj2->n6DOF;
    
    // Get object centers, velocities and angular velocities
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Relative position vectors from centers to contact point
    Eigen::Vector3d r1 = contact_point - center1;
    Eigen::Vector3d r2 = contact_point - center2;
    
    // Get linear velocities at centers
    Eigen::Vector3d v1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    Eigen::Vector3d v2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    
    // Get angular velocities
    Eigen::Vector3d omega1 = obj1->omega_I;
    Eigen::Vector3d omega2 = obj2->omega_I;
    
    // Calculate velocities at contact point
    Eigen::Vector3d v1_contact = v1 + omega1.cross(r1);
    Eigen::Vector3d v2_contact = v2 + omega2.cross(r2);
    
    // Relative velocity at contact point
    Eigen::Vector3d v_rel = v2_contact - v1_contact;
    
    // Normal component of relative velocity
    double v_rel_n = v_rel.dot(normal);
    
    // Tangential component of relative velocity
    Eigen::Vector3d v_rel_t = v_rel - v_rel_n * normal;
    double v_rel_t_mag = v_rel_t.norm();
    
    // Unit vector in tangential direction
    Eigen::Vector3d t_hat;
    if(v_rel_t_mag > 1.0e-10)
    {
        t_hat = v_rel_t / v_rel_t_mag;
    }
    else
    {
        // If tangential velocity is close to zero, use a default tangential direction
        // (perpendicular to normal)
        if(fabs(normal(0)) > 0.5)
            t_hat = Eigen::Vector3d(0.0, 1.0, 0.0);
        else
            t_hat = Eigen::Vector3d(1.0, 0.0, 0.0);
            
        t_hat = t_hat - normal.dot(t_hat) * normal;
        t_hat.normalize();
    }
    
    // Calculate normal force using linear spring-dashpot model
    double fn = spring_constant * overlap - damping_constant * v_rel_n;
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Calculate tangential (friction) force
    double ft = friction_coefficient * fn;
    if(v_rel_t_mag > 1.0e-10)
    {
        ft = min(ft, v_rel_t_mag); // Limit friction to prevent sticking
    }
    else
    {
        ft = 0.0;
    }
    
    // Total force vector
    force = fn * normal - ft * t_hat;
    
    // Calculate torque
    torque = r1.cross(force);
}

double sixdof_collision::calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2)
{
    // Get object centers
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Calculate distance between centers
    return (center2 - center1).norm();
} 