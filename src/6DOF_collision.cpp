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

#include"6DOF_collision.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"6DOF_collision_grid.h"
#include<math.h>
#include<iostream>

sixdof_collision::sixdof_collision(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF Collision Model startup..."<<endl;
    
    // Set default collision model
    contact_model = ContactForceModel::Linear;
    
    // Initialize common parameters from control file parameters
    // These should be added to the parameter file in a real implementation
    // For now, we set default values
    spring_constant_n = 1.0e6;            // Default normal stiffness [N/m]
    spring_constant_t = 0.5e6;            // Default tangential stiffness [N/m]
    damping_constant_n = 1.0e4;           // Default normal damping [N·s/m]
    damping_constant_t = 0.5e4;           // Default tangential damping [N·s/m]
    friction_coefficient = 0.3;           // Default friction coefficient
    restitution_coefficient = 0.7;        // Default restitution coefficient
    
    // Initialize rolling friction parameters
    rolling_friction_coefficient = 0.1;    // Default rolling friction coefficient
    rolling_stiffness = 0.5e6;             // Default rolling stiffness [N·m/rad]
    rolling_damping = 0.5e4;               // Default rolling damping [N·m·s/rad]
    rolling_torque_threshold = 1.0e-3;     // Default rolling torque threshold [N·m]
    
    // Initialize Hertzian contact parameters
    young_modulus = 1.0e7;                // Default Young's modulus [Pa]
    poisson_ratio = 0.3;                  // Default Poisson's ratio
    
    // Initialize DMT model parameters
    surface_energy = 0.05;                // Default surface energy [J/m²]
    dmt_cutoff_threshold = 0.1;           // Default cutoff threshold for DMT
    
    // Initialize JKR model parameters
    surface_energy_jkr = 0.1;             // Default surface energy for JKR [J/m²]
    jkr_cutoff_threshold = 0.2;           // Default cutoff threshold for JKR
    
    // Initialize sub-stepping parameters
    use_substeps = true;
    max_substeps = 10;
    
    // Create a new collision grid
    collision_grid = new sixdof_collision_grid(p, pgc);
}

sixdof_collision::~sixdof_collision()
{
    // Clean up collision grid
    if (collision_grid)
        delete collision_grid;
}

void sixdof_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Update contact history (remove pairs no longer in contact)
    update_contact_history(p);
    
    // Update the collision grid with current object positions
    collision_grid->update_grid(p, pgc, fb_obj);
    
    // Get potential collision pairs from the grid
    std::vector<std::pair<int, int>> potential_collisions = 
        collision_grid->find_potential_collisions(p, pgc, fb_obj);
    
    if(p->mpirank==0 && p->count%p->P12==0 && potential_collisions.size() > 0)
    {
        cout<<"6DOF Collision: Found "<<potential_collisions.size()<<" potential collision pairs"<<endl;
    }
    
    // Check each potential collision pair
    for(const auto& pair : potential_collisions)
    {
        int i = pair.first;
        int j = pair.second;
        
        // Variables to store collision information
        Eigen::Vector3d contact_point, normal;
        double overlap = 0.0;
        
        // Detect if collision occurred
        if(detect_collision(p, pgc, fb_obj[i], fb_obj[j], contact_point, normal, overlap))
        {
            // Calculate contact forces and torques
            Eigen::Vector3d force, torque, rolling_torque, twisting_torque;
            
            // Apply appropriate contact force model
            switch(contact_model)
            {
                case ContactForceModel::Linear:
                    calculate_linear_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                 contact_point, normal, overlap, 
                                                 force, torque);
                    break;
                    
                case ContactForceModel::Hertz:
                    calculate_hertz_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                contact_point, normal, overlap, 
                                                force, torque);
                    break;
                    
                case ContactForceModel::HertzMindlin:
                    calculate_hertz_mindlin_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                       contact_point, normal, overlap, 
                                                       force, torque);
                    break;
                    
                case ContactForceModel::DMT:
                    calculate_dmt_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                              contact_point, normal, overlap, 
                                              force, torque);
                    break;
                    
                case ContactForceModel::JKR:
                    calculate_jkr_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                              contact_point, normal, overlap, 
                                              force, torque);
                    break;
                    
                default:
                    // Default to linear model if unrecognized
                    calculate_linear_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                 contact_point, normal, overlap, 
                                                 force, torque);
                    break;
            }
            
            // Calculate rolling friction and twisting resistance
            calculate_rolling_friction_torque(p, pgc, fb_obj[i], fb_obj[j],
                                            contact_point, normal, overlap,
                                            rolling_torque);
                                            
            calculate_twisting_resistance(p, pgc, fb_obj[i], fb_obj[j],
                                        contact_point, normal, overlap,
                                        twisting_torque);
            
            // Add rolling friction and twisting resistance to total torque
            torque += rolling_torque + twisting_torque;
            
            // If using sub-stepping for collision resolution
            if(use_substeps && overlap > 0.01 * fb_obj[i]->radius) // Only use sub-stepping for significant overlaps
            {
                resolve_collision_with_substeps(p, pgc, fb_obj[i], fb_obj[j], 
                                             contact_point, normal, overlap, 
                                             force, torque);
            }
            else
            {
                // Add collision forces to object i (action)
                fb_obj[i]->Xext -= force(0);
                fb_obj[i]->Yext -= force(1);
                fb_obj[i]->Zext -= force(2);
                
                fb_obj[i]->Kext -= torque(0);
                fb_obj[i]->Mext -= torque(1);
                fb_obj[i]->Next -= torque(2);
                
                // Add collision forces to object j (reaction)
                fb_obj[j]->Xext += force(0);
                fb_obj[j]->Yext += force(1);
                fb_obj[j]->Zext += force(2);
                
                fb_obj[j]->Kext += torque(0);
                fb_obj[j]->Mext += torque(1);
                fb_obj[j]->Next += torque(2);
            }
            
            if(p->mpirank==0 && p->count%p->P12==0)
            {
                cout<<"6DOF Collision detected between objects "<<i<<" and "<<j<<endl;
                cout<<"  Model: ";
                switch(contact_model) {
                    case ContactForceModel::Linear: cout<<"Linear"; break;
                    case ContactForceModel::Hertz: cout<<"Hertz"; break;
                    case ContactForceModel::HertzMindlin: cout<<"Hertz-Mindlin"; break;
                    case ContactForceModel::DMT: cout<<"DMT"; break;
                    case ContactForceModel::JKR: cout<<"JKR"; break;
                    default: cout<<"Unknown"; break;
                }
                cout<<endl;
                cout<<"  Overlap: "<<overlap<<endl;
                cout<<"  Force: ["<<force(0)<<", "<<force(1)<<", "<<force(2)<<"]"<<endl;
            }
        }
    }
}

bool sixdof_collision::detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                      Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Using improved spherical collision detection as an initial check
    
    // Get object centers
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Calculate distance between centers
    Eigen::Vector3d center_diff = center2 - center1;
    double distance = center_diff.norm();
    
    // Calculate sum of bounding radii
    double sum_radii = obj1->radius + obj2->radius;
    
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
        contact_point = center1 + normal * (obj1->radius - 0.5 * overlap);
        
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
    
    // Get or create contact history for this pair
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end())
    {
        // New contact
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time;
    if(dt <= 0.0) dt = p->dt; // Use simulation dt if no meaningful history
    
    // Update contact history
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
    
    // Calculate normal force using linear spring-dashpot model
    // Fn = kn * delta - gamma_n * v_rel_n
    double fn = spring_constant_n * overlap - damping_constant_n * v_rel_n;
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Update tangential overlap based on relative velocity
    if(v_rel_t_mag > 1.0e-10)
    {
        // Increment tangential overlap
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap -= 
            normal * normal.dot(contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap);
        
        // Calculate tangential force based on spring and damping
        Eigen::Vector3d ft_vector = spring_constant_t * contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap 
                                  - damping_constant_t * v_rel_t;
        
        // Apply Coulomb's friction law
        double ft_mag = ft_vector.norm();
        double ft_max = friction_coefficient * fn;
        
        if(ft_mag > ft_max)
        {
            // Scale tangential force to the maximum allowed
            ft_vector *= (ft_max / ft_mag);
            
            // Update tangential overlap to match the maximum force
            contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap *= (ft_max / ft_mag);
        }
        
        // Total force vector
        force = fn * normal - ft_vector;
    }
    else
    {
        // No tangential motion, just apply normal force
        force = fn * normal;
    }
    
    // Calculate torque
    torque = r1.cross(force);
}

void sixdof_collision::calculate_hertz_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // Calculate effective radius and Young's modulus
    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = young_modulus / (2.0 * (1.0 - poisson_ratio * poisson_ratio));
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian model
    // Fn = (4/3) * k_hertz * delta^(3/2) - gamma_n * sqrt(delta) * v_rel_n
    double fn = (4.0/3.0) * k_hertz * pow(overlap, 1.5) - damping_constant_n * sqrt(overlap) * v_rel_n;
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Get or create contact history for this pair
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end())
    {
        // New contact
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time;
    if(dt <= 0.0) dt = p->dt; // Use simulation dt if no meaningful history
    
    // Update contact history
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
    
    // Update tangential overlap based on relative velocity
    if(v_rel_t_mag > 1.0e-10)
    {
        // Increment tangential overlap
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap -= 
            normal * normal.dot(contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap);
        
        // Calculate tangential force based on spring and damping
        // For Hertzian model, use sqrt(delta) scaling for tangential force
        Eigen::Vector3d ft_vector = spring_constant_t * sqrt(overlap) * contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap 
                                  - damping_constant_t * sqrt(overlap) * v_rel_t;
        
        // Apply Coulomb's friction law
        double ft_mag = ft_vector.norm();
        double ft_max = friction_coefficient * fn;
        
        if(ft_mag > ft_max)
        {
            // Scale tangential force to the maximum allowed
            ft_vector *= (ft_max / ft_mag);
            
            // Update tangential overlap to match the maximum force
            contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap *= (ft_max / ft_mag);
        }
        
        // Total force vector
        force = fn * normal - ft_vector;
    }
    else
    {
        // No tangential motion, just apply normal force
        force = fn * normal;
    }
    
    // Calculate torque
    torque = r1.cross(force);
}

void sixdof_collision::calculate_hertz_mindlin_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                         const Eigen::Vector3d &contact_point, 
                                                         const Eigen::Vector3d &normal, 
                                                         const double overlap,
                                                         Eigen::Vector3d &force, 
                                                         Eigen::Vector3d &torque)
{
    int id1 = obj1->n6DOF;
    int id2 = obj2->n6DOF;
    
    // Unique identifier for this contact pair
    pair<int, int> contact_pair = make_pair(min(id1, id2), max(id1, id2));
    
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
    
    // Get or create contact history for this pair
    auto it = contact_history.find(contact_pair);
    if(it == contact_history.end())
    {
        // New contact
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        contact_history[contact_pair] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[contact_pair].last_update_time;
    if(dt <= 0.0) dt = p->dt; // Use simulation dt if no meaningful history
    
    // Update contact history
    contact_history[contact_pair].in_contact = true;
    contact_history[contact_pair].last_update_time = p->simtime;
    
    // Calculate effective radius and effective Young's modulus
    double R_eff = calculate_effective_radius(obj1->radius, obj2->radius);
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    double fn = (4.0/3.0) * k_hertz * pow(overlap, 1.5) - damping_constant_n * v_rel_n;
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Calculate tangential spring constant (simplification of Mindlin theory)
    double G_eff = 0.5 * E_eff / (2.0 * (1.0 + poisson_ratio)); // Effective shear modulus
    double k_t = 8.0 * G_eff * sqrt(R_eff * overlap); // Tangential stiffness
    
    // Update tangential overlap (displacement)
    if(v_rel_t_mag > 1.0e-10)
    {
        Eigen::Vector3d t_hat = v_rel_t / v_rel_t_mag;
        
        // Increment tangential overlap based on relative velocity
        contact_history[contact_pair].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[contact_pair].tangential_overlap -= 
            normal * normal.dot(contact_history[contact_pair].tangential_overlap);
        
        // Calculate the tangential force based on tangential spring
        Eigen::Vector3d ft_vector = k_t * contact_history[contact_pair].tangential_overlap;
        double ft_mag = ft_vector.norm();
        
        // Apply Coulomb's friction law (capping the tangential force)
        double ft_max = friction_coefficient * fn;
        if(ft_mag > ft_max)
        {
            // Scale tangential overlap and force to the maximum allowed
            contact_history[contact_pair].tangential_overlap *= (ft_max / ft_mag);
            ft_vector *= (ft_max / ft_mag);
        }
        
        // Total force vector
        force = fn * normal - ft_vector;
    }
    else
    {
        // No tangential motion, just apply normal force
        force = fn * normal;
    }
    
    // Calculate torque
    torque = r1.cross(force);
}

void sixdof_collision::calculate_dmt_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // Calculate effective radius and Young's modulus
    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = young_modulus / (2.0 * (1.0 - poisson_ratio * poisson_ratio));
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate DMT force components
    // 1. Hertzian repulsive force
    double f_hertz = (4.0/3.0) * k_hertz * pow(overlap, 1.5);
    
    // 2. Van der Waals attractive force
    double f_vdw = 4.0 * M_PI * surface_energy * R_eff;
    
    // 3. Damping force
    double f_damp = damping_constant_n * sqrt(overlap) * v_rel_n;
    
    // Total normal force
    double fn = f_hertz - f_vdw - f_damp;
    
    // Apply cutoff threshold for DMT model
    if(overlap < dmt_cutoff_threshold)
    {
        fn = 0.0;
    }
    
    // Get or create contact history for this pair
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end())
    {
        // New contact
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time;
    if(dt <= 0.0) dt = p->dt; // Use simulation dt if no meaningful history
    
    // Update contact history
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
    
    // Update tangential overlap based on relative velocity
    if(v_rel_t_mag > 1.0e-10)
    {
        // Increment tangential overlap
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap -= 
            normal * normal.dot(contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap);
        
        // Calculate tangential force based on spring and damping
        // For DMT model, use sqrt(delta) scaling for tangential force
        Eigen::Vector3d ft_vector = spring_constant_t * sqrt(overlap) * contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap 
                                  - damping_constant_t * sqrt(overlap) * v_rel_t;
        
        // Apply Coulomb's friction law
        double ft_mag = ft_vector.norm();
        double ft_max = friction_coefficient * fn;
        
        if(ft_mag > ft_max)
        {
            // Scale tangential force to the maximum allowed
            ft_vector *= (ft_max / ft_mag);
            
            // Update tangential overlap to match the maximum force
            contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap *= (ft_max / ft_mag);
        }
        
        // Total force vector
        force = fn * normal - ft_vector;
    }
    else
    {
        // No tangential motion, just apply normal force
        force = fn * normal;
    }
    
    // Calculate torque
    torque = r1.cross(force);
}

double sixdof_collision::calculate_effective_young_modulus(double E1, double E2, double nu1, double nu2)
{
    // Calculate effective Young's modulus
    // 1/E* = (1-nu1²)/E1 + (1-nu2²)/E2
    double E_eff = 1.0 / ((1.0 - nu1*nu1)/E1 + (1.0 - nu2*nu2)/E2);
    return E_eff;
}

double sixdof_collision::calculate_effective_radius(double R1, double R2)
{
    // Calculate effective radius
    // 1/R* = 1/R1 + 1/R2
    double R_eff = 1.0 / (1.0/R1 + 1.0/R2);
    return R_eff;
}

double sixdof_collision::calculate_hertz_stiffness(double E_eff, double R_eff)
{
    // Return the effective stiffness for Hertzian contact
    // This is the coefficient in the formula: F = (4/3) * E* * sqrt(R*) * delta^(3/2)
    return E_eff * sqrt(R_eff);
}

void sixdof_collision::update_contact_history(lexer *p)
{
    // Remove contact history for pairs that are no longer in contact
    // and haven't been for a while
    double contact_timeout = 1.0; // seconds
    
    for(auto it = contact_history.begin(); it != contact_history.end();)
    {
        if(!it->second.in_contact && (p->simtime - it->second.last_update_time) > contact_timeout)
        {
            // Remove history for pairs no longer in contact
            it = contact_history.erase(it);
        }
        else
        {
            // Reset in_contact for this time step
            // It will be set to true if contact is detected
            it->second.in_contact = false;
            ++it;
        }
    }
}

double sixdof_collision::calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2)
{
    // Get object centers
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Calculate distance between centers
    return (center2 - center1).norm();
}

bool sixdof_collision::detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                               Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // First do a quick sphere-sphere check
    if(!detect_collision(p, pgc, obj1, obj2, contact_point, normal, overlap))
        return false;
        
    // If sphere-sphere check passes, do detailed triangle-triangle check
    double min_overlap = 1e10;
    bool collision_found = false;
    
    // Transform triangles to world coordinates
    for(int i=0; i<obj1->tricount; ++i)
    {
        // Get triangle vertices in world coordinates
        Eigen::Vector3d v1[3];
        for(int q=0; q<3; ++q)
        {
            Eigen::Vector3d local_point(obj1->tri_x0[i][q], obj1->tri_y0[i][q], obj1->tri_z0[i][q]);
            Eigen::Vector3d global_point;
            obj1->motionext_trans(p, pgc, local_point, global_point);
            v1[q] = global_point;
        }
        
        for(int j=0; j<obj2->tricount; ++j)
        {
            // Get triangle vertices in world coordinates
            Eigen::Vector3d v2[3];
            for(int q=0; q<3; ++q)
            {
                Eigen::Vector3d local_point(obj2->tri_x0[j][q], obj2->tri_y0[j][q], obj2->tri_z0[j][q]);
                Eigen::Vector3d global_point;
                obj2->motionext_trans(p, pgc, local_point, global_point);
                v2[q] = global_point;
            }
            
            // Check for triangle-triangle intersection
            Eigen::Vector3d contact, norm;
            double overlap_dist;
            
            if(triangle_triangle_intersection(v1, v2, contact, norm, overlap_dist))
            {
                if(overlap_dist < min_overlap)
                {
                    min_overlap = overlap_dist;
                    contact_point = contact;
                    normal = norm;
                    collision_found = true;
                }
            }
        }
    }
    
    if(collision_found)
    {
        overlap = min_overlap;
        return true;
    }
    
    return false;
}

bool sixdof_collision::triangle_triangle_intersection(const Eigen::Vector3d v1[3], const Eigen::Vector3d v2[3],
                                                    Eigen::Vector3d &contact, Eigen::Vector3d &normal, double &overlap)
{
    // Calculate triangle normals
    Eigen::Vector3d n1 = (v1[1] - v1[0]).cross(v1[2] - v1[0]).normalized();
    Eigen::Vector3d n2 = (v2[1] - v2[0]).cross(v2[2] - v2[0]).normalized();
    
    // Check if triangles are coplanar
    if(fabs(n1.dot(n2)) > 0.999)
        return false;
        
    // Calculate intersection line
    Eigen::Vector3d line_dir = n1.cross(n2).normalized();
    
    // Project vertices onto intersection line
    double t1[3], t2[3];
    for(int i=0; i<3; ++i)
    {
        t1[i] = v1[i].dot(line_dir);
        t2[i] = v2[i].dot(line_dir);
    }
    
    // Check for overlap in projections
    double min1 = std::min({t1[0], t1[1], t1[2]});
    double max1 = std::max({t1[0], t1[1], t1[2]});
    double min2 = std::min({t2[0], t2[1], t2[2]});
    double max2 = std::max({t2[0], t2[1], t2[2]});
    
    if(max1 < min2 || max2 < min1)
        return false;
        
    // Calculate overlap
    overlap = std::min(max1, max2) - std::max(min1, min2);
    
    // Calculate contact point (midpoint of overlap)
    double t_contact = 0.5 * (std::max(min1, min2) + std::min(max1, max2));
    contact = line_dir * t_contact;
    
    // Calculate normal (average of triangle normals)
    normal = (n1 + n2).normalized();
    
    return true;
}

void sixdof_collision::calculate_jkr_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // Calculate effective radius and Young's modulus
    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = young_modulus / (2.0 * (1.0 - poisson_ratio * poisson_ratio));
    
    // Calculate JKR parameters
    double gamma = surface_energy_jkr;  // Surface energy
    double a0 = pow(6.0 * M_PI * gamma * R_eff * R_eff / E_eff, 1.0/3.0);  // Contact radius at zero load
    
    // Calculate contact radius
    double a = sqrt(R_eff * overlap);  // Current contact radius
    
    // Calculate JKR force components
    // 1. Elastic force
    double f_elastic = (4.0/3.0) * E_eff * sqrt(R_eff) * pow(overlap, 1.5);
    
    // 2. Adhesive force
    double f_adhesive = -4.0 * M_PI * gamma * R_eff * (1.0 - pow(a0/a, 1.5));
    
    // 3. Damping force
    double f_damp = damping_constant_n * sqrt(overlap) * v_rel_n;
    
    // Total normal force
    double fn = f_elastic + f_adhesive - f_damp;
    
    // Apply cutoff threshold for JKR model
    if(overlap < jkr_cutoff_threshold)
    {
        fn = 0.0;
    }
    
    // Get or create contact history for this pair
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end())
    {
        // New contact
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time;
    if(dt <= 0.0) dt = p->dt; // Use simulation dt if no meaningful history
    
    // Update contact history
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
    
    // Update tangential overlap based on relative velocity
    if(v_rel_t_mag > 1.0e-10)
    {
        // Increment tangential overlap
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap -= 
            normal * normal.dot(contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap);
        
        // Calculate tangential force based on spring and damping
        // For JKR model, use sqrt(delta) scaling for tangential force
        Eigen::Vector3d ft_vector = spring_constant_t * sqrt(overlap) * contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap 
                                  - damping_constant_t * sqrt(overlap) * v_rel_t;
        
        // Apply Coulomb's friction law
        double ft_mag = ft_vector.norm();
        double ft_max = friction_coefficient * fn;
        
        if(ft_mag > ft_max)
        {
            // Scale tangential force to the maximum allowed
            ft_vector *= (ft_max / ft_mag);
            
            // Update tangential overlap to match the maximum force
            contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap *= (ft_max / ft_mag);
        }
        
        // Total force vector
        force = fn * normal - ft_vector;
    }
    else
    {
        // No tangential motion, just apply normal force
        force = fn * normal;
    }
    
    // Calculate torque
    torque = r1.cross(force);
}

void sixdof_collision::calculate_rolling_friction_torque(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                      const Eigen::Vector3d &contact_point,
                                                      const Eigen::Vector3d &normal,
                                                      const double overlap,
                                                      Eigen::Vector3d &rolling_torque)
{
    // Get object centers and angular velocities
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Relative position vectors from centers to contact point
    Eigen::Vector3d r1 = contact_point - center1;
    Eigen::Vector3d r2 = contact_point - center2;
    
    // Get angular velocities
    Eigen::Vector3d omega1 = obj1->omega_I;
    Eigen::Vector3d omega2 = obj2->omega_I;
    
    // Calculate relative angular velocity
    Eigen::Vector3d omega_rel = omega2 - omega1;
    
    // Project relative angular velocity onto the contact plane
    Eigen::Vector3d omega_rel_t = omega_rel - normal * normal.dot(omega_rel);
    
    // Calculate rolling velocity
    double rolling_velocity = omega_rel_t.norm();
    
    if(rolling_velocity > 1.0e-10)
    {
        // Calculate rolling direction
        Eigen::Vector3d rolling_direction = omega_rel_t / rolling_velocity;
        
        // Calculate rolling torque magnitude
        double rolling_torque_mag = rolling_friction_coefficient * overlap * 
                                  (rolling_stiffness * rolling_velocity + 
                                   rolling_damping * rolling_velocity * rolling_velocity);
        
        // Apply threshold
        if(rolling_torque_mag > rolling_torque_threshold)
        {
            rolling_torque_mag = rolling_torque_threshold;
        }
        
        // Calculate rolling torque vector
        rolling_torque = -rolling_torque_mag * rolling_direction;
    }
    else
    {
        rolling_torque.setZero();
    }
}

void sixdof_collision::calculate_twisting_resistance(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                  const Eigen::Vector3d &contact_point,
                                                  const Eigen::Vector3d &normal,
                                                  const double overlap,
                                                  Eigen::Vector3d &twisting_torque)
{
    // Get object centers and angular velocities
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    
    // Get angular velocities
    Eigen::Vector3d omega1 = obj1->omega_I;
    Eigen::Vector3d omega2 = obj2->omega_I;
    
    // Calculate relative angular velocity
    Eigen::Vector3d omega_rel = omega2 - omega1;
    
    // Project relative angular velocity onto the normal direction
    double twisting_velocity = normal.dot(omega_rel);
    
    if(std::abs(twisting_velocity) > 1.0e-10)
    {
        // Calculate twisting torque magnitude
        double twisting_torque_mag = rolling_friction_coefficient * overlap * 
                                   (rolling_stiffness * std::abs(twisting_velocity) + 
                                    rolling_damping * twisting_velocity * twisting_velocity);
        
        // Apply threshold
        if(twisting_torque_mag > rolling_torque_threshold)
        {
            twisting_torque_mag = rolling_torque_threshold;
        }
        
        // Calculate twisting torque vector
        twisting_torque = -twisting_torque_mag * normal * (twisting_velocity > 0 ? 1.0 : -1.0);
    }
    else
    {
        twisting_torque.setZero();
    }
} 