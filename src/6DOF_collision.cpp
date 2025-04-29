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
    
    // Set default collision model
    contact_model = ContactForceModel::Linear;
    
    // Initialize common parameters from control file parameters
    // These should be added to the parameter file in a real implementation
    // For now, we set default values
    spring_constant = 1.0e6;            // Default stiffness [N/m]
    damping_constant = 1.0e4;           // Default damping [N·s/m]
    friction_coefficient = 0.3;         // Default friction coefficient
    restitution_coefficient = 0.7;      // Default restitution coefficient
    
    // Initialize Hertzian contact parameters
    young_modulus = 1.0e7;              // Default Young's modulus [Pa]
    poisson_ratio = 0.3;                // Default Poisson's ratio
    
    // Initialize DMT model parameters
    surface_energy = 0.05;              // Default surface energy [J/m²]
    dmt_cutoff_threshold = 0.1;         // Default cutoff threshold for DMT
    
    // Allocate space for bounding radius for each 6DOF object
    bounding_radius.resize(p->X20);
    
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
    
    // Calculate effective mass
    // Following Lethe's approach for equivalent mass calculation
    double m1 = obj1->Mass_fb;
    double m2 = obj2->Mass_fb;
    double effective_mass = (m1 * m2) / (m1 + m2);
    
    // Calculate spring constant using effective mass and restitution coefficient
    // For now, we'll keep using the existing spring_constant parameter
    // but will add proper calculation in future updates
    
    // Calculate damping coefficient based on coefficient of restitution
    // Following Lethe's logarithmic approach
    double beta = 0.0;
    
    // Handle edge cases for restitution coefficient
    if(restitution_coefficient < 1.0e-6) 
    {
        // For near-zero restitution, use critical damping
        beta = 1.0;
    }
    else if(restitution_coefficient >= 1.0) 
    {
        // For perfect or super-elastic collisions, no damping
        beta = 0.0;
    }
    else
    {
        // Normal case: logarithmic model
        double log_coeff_restitution = log(restitution_coefficient);
        beta = log_coeff_restitution / sqrt(log_coeff_restitution * log_coeff_restitution + M_PI * M_PI);
    }
    
    // Calculate damping using the proper formula based on critical damping
    double critical_damping = 2.0 * sqrt(effective_mass * spring_constant);
    double damping = -beta * critical_damping;
    
    // Calculate normal force using linear spring-dashpot model with proper damping
    double fn = spring_constant * overlap;
    
    // Only apply damping when objects are approaching (v_rel_n < 0)
    if (v_rel_n < 0)
    {
        fn -= damping * v_rel_n;
    }
    
    fn = max(fn, 0.0); // Ensure normal force is repulsive
    
    // Calculate tangential (friction) force
    double ft = friction_coefficient * fn;
    if(v_rel_t_mag > 1.0e-10)
    {
        // Apply tangential force opposite to tangential velocity
        ft = min(ft, damping * v_rel_t_mag); // Limit friction force
    }
    else
    {
        ft = 0.0;
    }
    
    // Total force vector
    force = (fn * normal - ft * t_hat);
    
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
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    // Fn = (4/3) * E* * sqrt(R*) * delta^(3/2) - damping_coefficient * v_rel_n
    double fn = (4.0/3.0) * k_hertz * pow(overlap, 1.5) - damping_constant * v_rel_n;
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
    double R_eff = calculate_effective_radius(bounding_radius[id1], bounding_radius[id2]);
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    double fn = (4.0/3.0) * k_hertz * pow(overlap, 1.5) - damping_constant * v_rel_n;
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
    
    // Calculate effective radius and effective Young's modulus
    double R_eff = calculate_effective_radius(bounding_radius[id1], bounding_radius[id2]);
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    
    // Calculate Hertzian stiffness
    double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    
    // Calculate normal force using Hertzian contact model (non-linear spring with damping)
    double hertz_force = (4.0/3.0) * k_hertz * pow(overlap, 1.5);
    
    // Calculate adhesive force from DMT model
    double pull_off_force = 2.0 * M_PI * surface_energy * R_eff;
    
    // Apply contact force only if we're within the cutoff threshold
    // In the DMT model, adhesion forces act even at small separations
    double fn = hertz_force - pull_off_force;
    
    // Apply damping
    fn -= damping_constant * v_rel_n;
    
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