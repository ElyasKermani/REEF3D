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
#include<math.h>
#include<iostream>

sixdof_collision::sixdof_collision(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF Collision Model startup..."<<endl;
    
    // Set default collision model to Hertz-Mindlin
    contact_model = ContactForceModel::HertzMindlin;
    
    // Initialize base parameters (these are reference values for a 1m object)
    base_spring_constant = 1.0e7;            // Base stiffness [N/m]
    base_damping_constant = 2.0e4;           // Base damping [N·s/m]
    base_young_modulus = 2.0e7;              // Base Young's modulus [Pa]
    base_surface_energy = 0.1;               // Base surface energy [J/m²]
    
    // Initialize fixed parameters
    poisson_ratio = 0.3;                     // Poisson's ratio
    friction_coefficient = 0.3;              // Friction coefficient
    restitution_coefficient = 0.65;          // Restitution coefficient
    dmt_cutoff_threshold = 0.1;              // DMT cutoff threshold
    
    // Rolling resistance parameters
    rolling_friction_coefficient = 0.01;      // Rolling friction coefficient
    rolling_viscous_damping = 1.0e3;         // Rolling viscous damping [N·m·s]
    
    // Allocate space for bounding radius and contact history
    bounding_radius.resize(p->X20);
    contact_history.clear();
    
    // Set initial bounding radius for each object
    // In a real implementation, this would be calculated based on the object geometry
    for(int i=0; i<p->X20; ++i)
    {
        bounding_radius[i] = 0.12; // Default radius (should be calculated from geometry)
    }
}

sixdof_collision::~sixdof_collision()
{
}

void sixdof_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Update contact history (remove pairs no longer in contact)
    update_contact_history(p);

    // Check for collisions between all pairs of objects
    for(int i=0; i<p->X20-1; ++i)
    {
        for(int j=i+1; j<p->X20; ++j)
        {
            // Calculate size-dependent parameters for this pair
            calculate_size_dependent_parameters(p, fb_obj[i], fb_obj[j]);
            
            // Variables to store collision information
            Eigen::Vector3d contact_point, normal;
            double overlap = 0.0;
            
            // Detect if collision occurred
            if(detect_collision(p, pgc, fb_obj[i], fb_obj[j], contact_point, normal, overlap))
            {
                // Calculate forces and torques from collision
                Eigen::Vector3d force, torque;
                
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
                        
                    default:
                        // Default to linear model if unrecognized
                        calculate_linear_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                     contact_point, normal, overlap, 
                                                     force, torque);
                        break;
                }
                
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
                
                if(p->mpirank==0 && p->count%p->P12==0)
                {
                    cout<<"6DOF Collision detected between objects "<<i<<" and "<<j<<endl;
                    cout<<"  Model: ";
                    switch(contact_model) {
                        case ContactForceModel::Linear: cout<<"Linear"; break;
                        case ContactForceModel::Hertz: cout<<"Hertz"; break;
                        case ContactForceModel::HertzMindlin: cout<<"Hertz-Mindlin"; break;
                        case ContactForceModel::DMT: cout<<"DMT"; break;
                        default: cout<<"Unknown"; break;
                    }
                    cout<<endl;
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
    
    // Calculate effective mass for the collision
    double m1 = obj1->Mass_fb;
    double m2 = obj2->Mass_fb;
    double m_eff = (m1 * m2) / (m1 + m2);  // Reduced mass

    // Calculate effective stiffness based on Young's modulus and size
    double R1 = bounding_radius[obj1->n6DOF];
    double R2 = bounding_radius[obj2->n6DOF];
    double R_eff = (R1 * R2) / (R1 + R2);  // Effective radius
    
    // Calculate effective Young's modulus
    double E1 = scaled_young_modulus;
    double E2 = scaled_young_modulus;
    double nu1 = poisson_ratio;
    double nu2 = poisson_ratio;
    double E_eff = calculate_effective_young_modulus(E1, E2, nu1, nu2);
    
    // Calculate contact stiffness based on Hertz theory
    double k_hertz = (4.0/3.0) * E_eff * sqrt(R_eff);
    
    // Calculate normal force using linear spring-dashpot model
    double fn = k_hertz * overlap - scaled_damping_constant * v_rel_n;
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

    // Total force vector (no artificial scaling)
    force = fn * normal - ft * t_hat;

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
    
    // Unit vector in tangential direction
    Eigen::Vector3d t_hat;
    if(v_rel_t_mag > 1.0e-10)
    {
        t_hat = v_rel_t / v_rel_t_mag;
    }
    else
    {
        // If tangential velocity is close to zero, use a default tangential direction
        if(fabs(normal(0)) > 0.5)
            t_hat = Eigen::Vector3d(0.0, 1.0, 0.0);
        else
            t_hat = Eigen::Vector3d(1.0, 0.0, 0.0);
            
        t_hat = t_hat - normal.dot(t_hat) * normal;
        t_hat.normalize();
    }
    
    // Calculate effective radius and effective Young's modulus
    double R_eff = calculate_effective_radius(bounding_radius[id1], bounding_radius[id2]);
    double E_eff = calculate_effective_young_modulus(scaled_young_modulus, scaled_young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    // Fn = (4/3) * E* * sqrt(R*) * delta^(3/2) - damping_coefficient * v_rel_n
    double fn = (4.0/3.0) * k_hertz * pow(overlap, 1.5) - scaled_damping_constant * v_rel_n;
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Calculate tangential (friction) force - Coulomb friction
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

void sixdof_collision::calculate_hertz_mindlin_contact_force(lexer *p, ghostcell *pgc,
                                                           sixdof_obj *obj1, sixdof_obj *obj2,
                                                           const Eigen::Vector3d &contact_point,
                                                           const Eigen::Vector3d &normal,
                                                           const double overlap,
                                                           Eigen::Vector3d &force,
                                                           Eigen::Vector3d &torque)
{
    // Get object IDs for contact history
    int id1 = obj1->n6DOF;
    int id2 = obj2->n6DOF;
    auto contact_pair = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    
    // Get or create contact history
    auto &history = contact_history[contact_pair];
    if(!history.in_contact)
    {
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
    }
    
    // Calculate Hertz stiffness
    double k_hertz = calculate_hertz_stiffness(effective_young_modulus, effective_radius);
    
    // Get relative velocity at contact point
    Eigen::Vector3d vel1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    vel1 += obj1->omega_I.cross(contact_point - obj1->c_);
    Eigen::Vector3d vel2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    vel2 += obj2->omega_I.cross(contact_point - obj2->c_);
    Eigen::Vector3d rel_vel = vel2 - vel1;
    
    // Normal component of relative velocity
    double v_n = rel_vel.dot(normal);
    
    // Calculate normal force using Hertz-Mindlin model with force limiting
    double F_n = k_hertz * pow(overlap, 1.5);  // Hertz normal force
    double F_n_damping = scaled_damping_constant * v_n;  // Damping force
    
    // Limit maximum force based on material strength
    double max_force = M_PI * effective_radius * effective_radius * scaled_young_modulus;
    F_n = std::min(F_n + F_n_damping, max_force);
    
    // Calculate tangential force
    Eigen::Vector3d v_t = rel_vel - v_n * normal;
    double dt = p->simtime - history.last_update_time;
    history.last_update_time = p->simtime;
    
    // Update tangential overlap
    if(dt > 0.0)
    {
        history.tangential_overlap += v_t * dt;
        
        // Project tangential overlap to be perpendicular to normal
        history.tangential_overlap -= normal * history.tangential_overlap.dot(normal);
    }
    
    // Calculate tangential stiffness (Mindlin)
    double G = effective_young_modulus / (2.0 * (1.0 + poisson_ratio));
    double k_t = 8.0 * G * sqrt(effective_radius * overlap);
    
    // Calculate tangential force
    Eigen::Vector3d F_t = -k_t * history.tangential_overlap;
    double F_t_mag = F_t.norm();
    
    // Apply Coulomb friction limit
    double F_t_max = friction_coefficient * F_n;
    if(F_t_mag > F_t_max)
    {
        F_t *= F_t_max / F_t_mag;
        // Reset tangential overlap to match limited force
        history.tangential_overlap = -F_t / k_t;
    }
    
    // Calculate rolling resistance torque
    Eigen::Vector3d omega_rel = obj2->omega_I - obj1->omega_I;
    Eigen::Vector3d rolling_torque;
    if(omega_rel.norm() > 1e-10)
    {
        // Viscous rolling resistance
        rolling_torque = -scaled_rolling_damping * omega_rel;
        
        // Limit rolling torque based on normal force
        double max_rolling_torque = rolling_friction_coefficient * F_n * effective_radius;
        if(rolling_torque.norm() > max_rolling_torque)
        {
            rolling_torque *= max_rolling_torque / rolling_torque.norm();
        }
    }
    else
    {
        rolling_torque.setZero();
    }
    
    // Combine forces
    force = F_n * normal + F_t;
    
    // Calculate torque from force and rolling resistance
    Eigen::Vector3d r1 = contact_point - obj1->c_;
    torque = r1.cross(force) + rolling_torque;
    
    // Log forces if debugging
    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"  Normal force: "<<F_n<<endl;
        cout<<"  Tangential force: "<<F_t.norm()<<endl;
        cout<<"  Rolling torque: "<<rolling_torque.norm()<<endl;
        cout<<"  Effective radius: "<<effective_radius<<endl;
        cout<<"  Effective Young's modulus: "<<effective_young_modulus<<endl;
    }
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
    
    // Calculate effective radius and effective Young's modulus
    double R_eff = calculate_effective_radius(bounding_radius[id1], bounding_radius[id2]);
    double E_eff = calculate_effective_young_modulus(scaled_young_modulus, scaled_young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    double hertz_force = (4.0/3.0) * k_hertz * pow(overlap, 1.5);
    
    // Calculate adhesive force from DMT model
    double pull_off_force = 2.0 * M_PI * scaled_surface_energy * R_eff;
    
    // Apply contact force only if we're within the cutoff threshold
    // In the DMT model, adhesion forces act even at small separations
    double fn = hertz_force - pull_off_force;
    
    // Apply damping
    fn -= scaled_damping_constant * v_rel_n;
    
    // Handle friction
    Eigen::Vector3d ft_vector = Eigen::Vector3d::Zero();
    if(v_rel_t_mag > 1.0e-10)
    {
        Eigen::Vector3d t_hat = v_rel_t / v_rel_t_mag;
        
        // Calculate the maximum tangential force (Coulomb friction)
        double ft_max = friction_coefficient * fabs(fn);
        
        // Apply tangential force
        ft_vector = -ft_max * t_hat;
    }
    
    // Total force vector
    force = fn * normal + ft_vector;
    
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

void sixdof_collision::calculate_size_dependent_parameters(lexer *p, sixdof_obj *obj1, sixdof_obj *obj2)
{
    // Get characteristic sizes (using bounding radius as approximation)
    double size1 = bounding_radius[obj1->n6DOF] * 2.0;  // Diameter
    double size2 = bounding_radius[obj2->n6DOF] * 2.0;  // Diameter
    
    // Calculate average size ratio compared to base 1m object
    double size_ratio = (size1 + size2) / 2.0;
    
    // Scale spring constant linearly with size
    scaled_spring_constant = base_spring_constant * size_ratio;
    
    // Scale damping with square root of size (to maintain critical damping ratio)
    scaled_damping_constant = base_damping_constant * sqrt(size_ratio);
    
    // Scale Young's modulus linearly with size
    scaled_young_modulus = base_young_modulus * size_ratio;
    
    // Scale surface energy with square of size (surface area scaling)
    scaled_surface_energy = base_surface_energy * size_ratio * size_ratio;
    
    // Calculate effective material properties
    double E1 = scaled_young_modulus;
    double E2 = scaled_young_modulus;
    double nu1 = poisson_ratio;
    double nu2 = poisson_ratio;
    
    effective_young_modulus = calculate_effective_young_modulus(E1, E2, nu1, nu2);
    effective_radius = calculate_effective_radius(bounding_radius[obj1->n6DOF], 
                                                bounding_radius[obj2->n6DOF]);
    
    // Scale rolling resistance parameters
    scaled_rolling_damping = rolling_viscous_damping * size_ratio * size_ratio;
    
    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"6DOF Collision parameters scaled for size: "<<size_ratio<<"m"<<endl;
        cout<<"  Spring constant: "<<scaled_spring_constant<<" N/m"<<endl;
        cout<<"  Damping constant: "<<scaled_damping_constant<<" N·s/m"<<endl;
        cout<<"  Young's modulus: "<<scaled_young_modulus<<" Pa"<<endl;
        cout<<"  Surface energy: "<<scaled_surface_energy<<" J/m²"<<endl;
    }
} 