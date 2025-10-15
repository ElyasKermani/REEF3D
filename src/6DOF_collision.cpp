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
    
    // Set maximum number of objects and initialize storage
    max_objects = p->X20;
    collision_forces.resize(max_objects);
    collision_torques.resize(max_objects);
    object_aabbs.resize(max_objects);  // Allocate AABB storage
    clear_collision_forces();
    
    // Set default collision model
    contact_model = ContactForceModel::PhasicFlowNonLinearNonLimited;
    
    // Initialize common parameters from control file parameters
    // These should be added to the parameter file in a real implementation
    // For now, we set default values
    spring_constant_n = 1.0e6;            // Default normal stiffness [N/m]
    spring_constant_t = 0.5e6;            // Default tangential stiffness [N/m]
    damping_constant_n = 1.0e4;           // Default normal damping [N·s/m]
    damping_constant_t = 0.5e4;           // Default tangential damping [N·s/m]
    friction_coefficient = 0.3;           // Default friction coefficient
    restitution_coefficient = 0.8;        // Default restitution coefficient
    
    // Initialize rolling friction parameters
    rolling_friction_coefficient = 0.2;    // Default rolling friction coefficient
    rolling_stiffness = 0.5e6;             // Default rolling stiffness [N·m/rad]
    rolling_damping = 0.5e4;               // Default rolling damping [N·m·s/rad]
    rolling_torque_threshold = 1.0e-3;     // Default rolling torque threshold [N·m]
    
    // Initialize Hertzian contact parameters
    young_modulus = 1.0e6;                // Default Young's modulus [Pa]
    poisson_ratio = 0.25;                  // Default Poisson's ratio

    // Initialize Linear model restitution-based damping options (disabled by default)
    linear_use_restitution = true;
    linear_en = restitution_coefficient;   // Default to same as global restitution
    linear_et = -1.0;                      // Use linear_en if not provided
    
    // Initialize DMT model parameters
    surface_energy = 0.05;                // Default surface energy [J/m²]
    dmt_cutoff_threshold = 0.1;           // Default cutoff threshold for DMT
    
    // Initialize JKR model parameters
    surface_energy_jkr = 0.1;             // Default surface energy for JKR [J/m²]
    jkr_cutoff_threshold = 0.2;           // Default cutoff threshold for JKR
    
    // NEW: Initialize PacIFiC-style enhanced parameters
    pacific_Es = 1.0e7;                   // Effective Young's modulus [Pa]
    pacific_en = 0.7;                     // Normal restitution coefficient
    pacific_Gs = 4.0e6;                   // Effective shear modulus [Pa] (≈ 4*Es for most materials)
    pacific_muc = 0.3;                    // Coulomb friction coefficient
    pacific_kr = 0.1;                     // Rolling resistance coefficient [s/m]
    pacific_beta = log(pacific_en) / sqrt(M_PI * M_PI + log(pacific_en) * log(pacific_en));
    pacific_m2sqrt56 = -2.0 * sqrt(5.0 / 6.0);
    
    // Initialize Hooke model parameters (PacIFiC style)
    hooke_kn = 1.0e6;                     // Normal stiffness coefficient [N/m]
    hooke_kt = 0.5e6;                     // Tangential stiffness coefficient [N/m]
    hooke_en = 0.7;                       // Normal restitution coefficient
    hooke_et = -1.0;                      // Tangential restitution coefficient (-1 for auto-calculation)
    hooke_muc = 0.3;                      // Coulomb friction coefficient
    hooke_kr = 0.1;                       // Rolling resistance coefficient [s/m]
    
    // Initialize sub-stepping parameters
    use_substeps = true;
    max_substeps = 10;
    
    // HYBRID: Initialize adaptive collision detection parameters
    use_adaptive_collision = true;           // Enable by default
    collision_simple_threshold = 50;         // < 50 triangles: use sphere-sphere (very fast)
    collision_moderate_threshold = 500;      // < 500 triangles: use triangle-triangle with BVH
    // >= 500 triangles: use triangle-triangle with BVH (still fast due to BVH)
    
    if(p->mpirank==0)
    {
        cout<<"6DOF Collision: Adaptive algorithm selection enabled"<<endl;
        cout<<"  Simple threshold (sphere-sphere): < "<<collision_simple_threshold<<" triangles"<<endl;
        cout<<"  Moderate threshold (triangle-mesh): < "<<collision_moderate_threshold<<" triangles"<<endl;
    }
    
    // Create a new collision grid
    collision_grid = new sixdof_collision_grid(p, pgc);
}

sixdof_collision::~sixdof_collision()
{
    // Clean up collision grid
    if (collision_grid)
        delete collision_grid;
}

void sixdof_collision::update_aabbs(vector<sixdof_obj*> &fb_obj)
{
    // Update AABB for each object based on current position and radius
    for(int i = 0; i < max_objects; ++i)
    {
        object_aabbs[i].update_from_sphere(fb_obj[i]->c_, fb_obj[i]->radius);
    }
}

void sixdof_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Clear collision forces for all objects on all processors
    clear_collision_forces();
    
    // Only rank 0 calculates collision forces
    if(p->mpirank == 0)
    {
        // Update AABBs for fast pre-filtering
        update_aabbs(fb_obj);
        
        // Update contact history (remove pairs no longer in contact)
        update_contact_history(p);
        
        // Update the collision grid with current object positions
        collision_grid->update_grid(p, pgc, fb_obj);
        
        // Get potential collision pairs from the grid
        std::vector<std::pair<int, int>> potential_collisions = 
            collision_grid->find_potential_collisions(p, pgc, fb_obj);
        
        if(p->count%p->P12==0 && potential_collisions.size() > 0)
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
            
            // AABB pre-check for fast rejection (avoids expensive collision tests)
            if(!object_aabbs[i].overlaps(object_aabbs[j]))
            {
                continue;  // AABBs don't overlap, skip this pair
            }
            
            // Detect collision using adaptive algorithm selection
            bool collision_detected = false;
            
            if(use_adaptive_collision)
            {
                // HYBRID: Automatically choose best algorithm based on complexity
                collision_detected = detect_collision_adaptive(p, pgc, fb_obj[i], fb_obj[j], 
                                                              contact_point, normal, overlap);
            }
            else
            {
                // Fallback: Use original sphere-sphere method
                collision_detected = detect_collision(p, pgc, fb_obj[i], fb_obj[j], 
                                                     contact_point, normal, overlap);
            }
            
            if(collision_detected)
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
                    case ContactForceModel::PhasicFlowLinearLimited:
                        calculate_phasicflow_linear_limited(p, pgc, fb_obj[i], fb_obj[j],
                                                           contact_point, normal, overlap,
                                                           force, torque);
                        break;
                    case ContactForceModel::PhasicFlowLinearNonLimited:
                        calculate_phasicflow_linear_nonlimited(p, pgc, fb_obj[i], fb_obj[j],
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

                    case ContactForceModel::PhasicFlowNonLinearLimited:
                        calculate_phasicflow_nonlinear_limited(p, pgc, fb_obj[i], fb_obj[j],
                                                               contact_point, normal, overlap,
                                                               force, torque);
                        break;

                    case ContactForceModel::PhasicFlowNonLinearNonLimited:
                        calculate_phasicflow_nonlinear_nonlimited(p, pgc, fb_obj[i], fb_obj[j],
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
                        
                    case ContactForceModel::PacIFiCHertz:
                        calculate_pacific_hertz_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
                                                           contact_point, normal, overlap, 
                                                           force, torque);
                        break;
                        
                    case ContactForceModel::PacIFiCHooke:
                        calculate_pacific_hooke_contact_force(p, pgc, fb_obj[i], fb_obj[j], 
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
                
                // Store collision forces and torques (instead of directly applying to Xext, etc.)
                // Object i gets negative force (action)
                collision_forces[i] -= force;
                collision_torques[i] -= torque;
                
                // Object j gets positive force (reaction)
                collision_forces[j] += force;
                collision_torques[j] += torque;
                
                if(p->count%p->P12==0)
                {
                    cout<<"6DOF Collision detected between objects "<<i<<" and "<<j<<endl;
                    cout<<"  Model: ";
                    switch(contact_model) {
                        case ContactForceModel::Linear: cout<<"Linear"; break;
                        case ContactForceModel::PhasicFlowLinearLimited: cout<<"PhasicFlow-Linear-Limited"; break;
                        case ContactForceModel::PhasicFlowLinearNonLimited: cout<<"PhasicFlow-Linear-NonLimited"; break;
                        case ContactForceModel::Hertz: cout<<"Hertz"; break;
                        case ContactForceModel::HertzMindlin: cout<<"Hertz-Mindlin"; break;
                        case ContactForceModel::DMT: cout<<"DMT"; break;
                        case ContactForceModel::JKR: cout<<"JKR"; break;
                        case ContactForceModel::PacIFiCHertz: cout<<"PacIFiC-Hertz"; break;
                        case ContactForceModel::PacIFiCHooke: cout<<"PacIFiC-Hooke"; break;
                        case ContactForceModel::PhasicFlowNonLinearLimited: cout<<"PhasicFlow-NonLinear-Limited"; break;
                        case ContactForceModel::PhasicFlowNonLinearNonLimited: cout<<"PhasicFlow-NonLinear-NonLimited"; break;
                        default: cout<<"Unknown"; break;
                    }
                    cout<<endl;
                    cout<<"  Overlap: "<<overlap<<endl;
                    cout<<"  Force: ["<<force(0)<<", "<<force(1)<<", "<<force(2)<<"]"<<endl;
                }
            }
        }
    }
    
    // Broadcast collision forces and torques from rank 0 to all processors
    broadcast_collision_forces(p, pgc);
    
    // Apply collision forces to external force variables on ALL processors
    for(int i = 0; i < max_objects; ++i)
    {
        fb_obj[i]->Xext += collision_forces[i](0);
        fb_obj[i]->Yext += collision_forces[i](1);
        fb_obj[i]->Zext += collision_forces[i](2);
        
        fb_obj[i]->Kext += collision_torques[i](0);
        fb_obj[i]->Mext += collision_torques[i](1);
        fb_obj[i]->Next += collision_torques[i](2);
    }
    
    // Debug output to verify MPI communication (only print occasionally)
    if(p->count%p->P12==0)
    {
        if(p->mpirank == 0)
        {
            cout<<"6DOF Collision: Applied forces to all objects on rank 0:"<<endl;
            for(int i = 0; i < max_objects; ++i)
            {
                if(collision_forces[i].norm() > 1.0e-10 || collision_torques[i].norm() > 1.0e-10)
                {
                    cout<<"  Object "<<i<<": Force=["<<collision_forces[i](0)<<", "<<collision_forces[i](1)<<", "<<collision_forces[i](2)<<"]"<<endl;
                    cout<<"            Torque=["<<collision_torques[i](0)<<", "<<collision_torques[i](1)<<", "<<collision_torques[i](2)<<"]"<<endl;
                }
            }
        }
        else
        {
            // Verify that non-rank-0 processors received the forces
            cout<<"6DOF Collision: Rank "<<p->mpirank<<" received forces for "<<max_objects<<" objects"<<endl;
        }
        
        // Uncomment the next line for debugging MPI synchronization (adds computational overhead)
        // verify_collision_forces_synchronization(p, pgc);
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
    // Optionally compute damping from restitution: c = 2*zeta*sqrt(k*meff)
    double c_n = damping_constant_n;
    if (linear_use_restitution) {
        const double en_lin = std::min(std::max(linear_en, 1.0e-6), 0.999999);
        const double zeta = -std::log(en_lin) / std::sqrt(M_PI*M_PI + std::log(en_lin)*std::log(en_lin));
        const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
        c_n = 2.0 * zeta * std::sqrt(spring_constant_n * meff);
    }
    // Fn = kn * delta - c_n * v_rel_n
    double fn = spring_constant_n * overlap - c_n * v_rel_n;
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
        double c_t = damping_constant_t;
        if (linear_use_restitution) {
            double et_lin;
            if (linear_et == -1.0) {
                const double en_lin = std::min(std::max(linear_en, 1.0e-6), 0.999999);
                et_lin = en_lin;
            } else {
                et_lin = std::min(std::max(linear_et, 1.0e-6), 0.999999);
            }
            const double zeta_t = -std::log(et_lin) / std::sqrt(M_PI*M_PI + std::log(et_lin)*std::log(et_lin));
            const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
            c_t = 2.0 * zeta_t * std::sqrt(spring_constant_t * meff);
        }
        Eigen::Vector3d ft_vector = spring_constant_t * contact_history[std::make_pair(std::min(id1, id2), std::max(id1, id2))].tangential_overlap 
                                  - c_t * v_rel_t;
        
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

// phasicFlow linear model with history limiting
void sixdof_collision::calculate_phasicflow_linear_limited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                      const Eigen::Vector3d &contact_point,
                                                      const Eigen::Vector3d &normal,
                                                      const double overlap,
                                                      Eigen::Vector3d &force,
                                                      Eigen::Vector3d &torque)
{
    // Use the same structure as calculate_linear_contact_force but enforce
    // phasicFlow's history rescaling behavior when Coulomb limit is hit and
    // compute damping constants directly from restitution if linear_use_restitution is true.

    int id1 = obj1->n6DOF;
    int id2 = obj2->n6DOF;

    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;
    Eigen::Vector3d r1 = contact_point - center1;
    Eigen::Vector3d r2 = contact_point - center2;

    Eigen::Vector3d v1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    Eigen::Vector3d v2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    Eigen::Vector3d omega1 = obj1->omega_I;
    Eigen::Vector3d omega2 = obj2->omega_I;

    Eigen::Vector3d v1_contact = v1 + omega1.cross(r1);
    Eigen::Vector3d v2_contact = v2 + omega2.cross(r2);
    Eigen::Vector3d v_rel = v2_contact - v1_contact;

    double v_rel_n = v_rel.dot(normal);
    Eigen::Vector3d v_rel_t = v_rel - v_rel_n * normal;

    auto key = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    auto it = contact_history.find(key);
    if(it == contact_history.end())
    {
        ContactHistory history; history.tangential_overlap.setZero(); history.in_contact = true; history.last_update_time = p->simtime;
        contact_history[key] = history;
    }
    double dt = p->simtime - contact_history[key].last_update_time; if(dt <= 0.0) dt = p->dt;
    contact_history[key].in_contact = true; contact_history[key].last_update_time = p->simtime;

    // phasicFlow computes damping from restitution inputs; we mirror linear_use_restitution option
    double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
    double c_n = damping_constant_n;
    double c_t = damping_constant_t;
    if (linear_use_restitution)
    {
        const double en_lin = std::min(std::max(linear_en, 1.0e-6), 0.999999);
        const double zeta_n = -std::log(en_lin) / std::sqrt(M_PI*M_PI + std::log(en_lin)*std::log(en_lin));
        c_n = 2.0 * zeta_n * std::sqrt(spring_constant_n * meff);
        double et_lin = (linear_et == -1.0) ? en_lin : std::min(std::max(linear_et, 1.0e-6), 0.999999);
        const double zeta_t = -std::log(et_lin) / std::sqrt(M_PI*M_PI + std::log(et_lin)*std::log(et_lin));
        c_t = 2.0 * zeta_t * std::sqrt(spring_constant_t * meff);
    }

    double fn = spring_constant_n * overlap - c_n * v_rel_n;
    fn = std::max(fn, 0.0);

    // Update tangential history and force
    if(v_rel_t.norm() > 1.0e-10)
    {
        auto &overlap_t = contact_history[key].tangential_overlap;
        overlap_t += v_rel_t * dt;
        overlap_t -= normal * normal.dot(overlap_t);

        Eigen::Vector3d Ft = -spring_constant_t * overlap_t - c_t * v_rel_t;
        double ft = Ft.norm();
        double ft_max = friction_coefficient * fn;
        if(ft > ft_max)
        {
            if(overlap_t.norm() > 0.0)
            {
                Ft *= (ft_max/ft);
                // limited: rescale history to match capped force
                overlap_t = -Ft / spring_constant_t;
            }
            else
            {
                Ft.setZero();
            }
        }
        force = fn * normal + Ft; // Ft already negative along tangential direction
    }
    else
    {
        force = fn * normal;
    }
    torque = r1.cross(force);
}

// phasicFlow linear model without history limiting
void sixdof_collision::calculate_phasicflow_linear_nonlimited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                         const Eigen::Vector3d &contact_point,
                                                         const Eigen::Vector3d &normal,
                                                         const double overlap,
                                                         Eigen::Vector3d &force,
                                                         Eigen::Vector3d &torque)
{
    // Same as limited but do not rescale history when Coulomb cap is hit (force capped only)
    int id1 = obj1->n6DOF; int id2 = obj2->n6DOF;
    Eigen::Vector3d center1 = obj1->c_; Eigen::Vector3d center2 = obj2->c_;
    Eigen::Vector3d r1 = contact_point - center1; Eigen::Vector3d r2 = contact_point - center2;
    Eigen::Vector3d v1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    Eigen::Vector3d v2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    Eigen::Vector3d omega1 = obj1->omega_I; Eigen::Vector3d omega2 = obj2->omega_I;
    Eigen::Vector3d v1_contact = v1 + omega1.cross(r1); Eigen::Vector3d v2_contact = v2 + omega2.cross(r2);
    Eigen::Vector3d v_rel = v2_contact - v1_contact; double v_rel_n = v_rel.dot(normal);
    Eigen::Vector3d v_rel_t = v_rel - v_rel_n * normal;

    auto key = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    if(contact_history.find(key) == contact_history.end()) { ContactHistory h; h.tangential_overlap.setZero(); h.in_contact=true; h.last_update_time=p->simtime; contact_history[key]=h; }
    double dt = p->simtime - contact_history[key].last_update_time; if(dt<=0.0) dt=p->dt;
    contact_history[key].in_contact=true; contact_history[key].last_update_time=p->simtime;

    double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
    double c_n = damping_constant_n; double c_t = damping_constant_t;
    if (linear_use_restitution)
    {
        const double en_lin = std::min(std::max(linear_en, 1.0e-6), 0.999999);
        const double zeta_n = -std::log(en_lin) / std::sqrt(M_PI*M_PI + std::log(en_lin)*std::log(en_lin));
        c_n = 2.0 * zeta_n * std::sqrt(spring_constant_n * meff);
        double et_lin = (linear_et == -1.0) ? en_lin : std::min(std::max(linear_et, 1.0e-6), 0.999999);
        const double zeta_t = -std::log(et_lin) / std::sqrt(M_PI*M_PI + std::log(et_lin)*std::log(et_lin));
        c_t = 2.0 * zeta_t * std::sqrt(spring_constant_t * meff);
    }

    double fn = spring_constant_n * overlap - c_n * v_rel_n; fn = std::max(fn, 0.0);
    if(v_rel_t.norm() > 1.0e-10)
    {
        auto &overlap_t = contact_history[key].tangential_overlap;
        overlap_t += v_rel_t * dt; overlap_t -= normal * normal.dot(overlap_t);
        Eigen::Vector3d Ft = -spring_constant_t * overlap_t - c_t * v_rel_t;
        double ft = Ft.norm(); double ft_max = friction_coefficient * fn;
        if(ft > ft_max)
        {
            if(overlap_t.norm() > 0.0) { Ft *= (ft_max/ft); }
            else { Ft.setZero(); }
        }
        force = fn * normal + Ft;
    }
    else { force = fn * normal; }
    torque = r1.cross(force);
}

// phasicFlow non-linear Hertz–Mindlin with no tangential damping, limited variant
void sixdof_collision::calculate_phasicflow_nonlinear_limited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                         const Eigen::Vector3d &contact_point,
                                                         const Eigen::Vector3d &normal,
                                                         const double overlap,
                                                         Eigen::Vector3d &force,
                                                         Eigen::Vector3d &torque)
{
    int id1 = obj1->n6DOF; int id2 = obj2->n6DOF;
    Eigen::Vector3d center1 = obj1->c_; Eigen::Vector3d center2 = obj2->c_;
    Eigen::Vector3d r1 = contact_point - center1; Eigen::Vector3d r2 = contact_point - center2;
    Eigen::Vector3d v1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    Eigen::Vector3d v2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    Eigen::Vector3d omega1 = obj1->omega_I; Eigen::Vector3d omega2 = obj2->omega_I;
    Eigen::Vector3d v1_contact = v1 + omega1.cross(r1); Eigen::Vector3d v2_contact = v2 + omega2.cross(r2);
    Eigen::Vector3d v_rel = v2_contact - v1_contact; double v_rel_n = v_rel.dot(normal);
    Eigen::Vector3d v_rel_t = v_rel - v_rel_n * normal;

    auto key = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    if(contact_history.find(key) == contact_history.end()) { ContactHistory h; h.tangential_overlap.setZero(); h.in_contact=true; h.last_update_time=p->simtime; contact_history[key]=h; }
    double dt = p->simtime - contact_history[key].last_update_time; if(dt<=0.0) dt=p->dt;
    contact_history[key].in_contact=true; contact_history[key].last_update_time=p->simtime;

    // Effective properties (mirror phasicFlow symbols: Yeff, Geff), use REEF3D young_modulus/poisson_ratio
    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    double G = young_modulus / (1.0 + poisson_ratio);
    //double G = young_modulus / (2.0 * (1.0 + poisson_ratio));
    // phasicFlow uses Geff directly; we approximate with same material => G

    // Normal elastic and damping per phasicFlow nonLinearCF: Fn = -(4/3)Y*√R δ^{3/2} - sqrt(meff*K_hz)*etha_n*δ^{1/4} v_n
    double mi = obj1->Mass_fb; double mj = obj2->Mass_fb; double meff = (mi*mj)/(mi+mj);
    double K_hertz = (4.0/3.0) * E_eff * std::sqrt(R_eff);
    // phasicFlow: ethan = -2.2664*ln(en)/sqrt(ln^2+π^2)
    const double en = std::min(std::max(restitution_coefficient, 1.0e-6), 0.999999);
    double ethan = -2.2664 * std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));
    Eigen::Vector3d Fn = (-(4.0/3.0) * E_eff * std::sqrt(R_eff) * std::pow(overlap, 1.5)
                          - std::sqrt(meff * K_hertz) * ethan * std::pow(overlap, 0.25) * v_rel_n) * normal;

    // Tangential elastic only (no tangential damping), history-based
    auto &overlap_t = contact_history[key].tangential_overlap;
    if(v_rel_t.norm() > 1.0e-10)
    {
        overlap_t += v_rel_t * dt; overlap_t -= normal * normal.dot(overlap_t);
        double kt = 8.0 * G * std::sqrt(R_eff * overlap);
        Eigen::Vector3d Ft = - kt * overlap_t;
        double ft = Ft.norm(); double ft_max = friction_coefficient * Fn.norm();
        if(ft > ft_max)
        {
            if(overlap_t.norm() > 0.0)
            {
                Ft *= (ft_max/ft);
                overlap_t = -Ft / kt; // limited: rescale
            }
            else { Ft.setZero(); }
        }
        force = Fn + Ft;
    }
    else { force = Fn; }
    torque = r1.cross(force);
}

// phasicFlow non-linear without history limiting
void sixdof_collision::calculate_phasicflow_nonlinear_nonlimited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                            const Eigen::Vector3d &contact_point,
                                                            const Eigen::Vector3d &normal,
                                                            const double overlap,
                                                            Eigen::Vector3d &force,
                                                            Eigen::Vector3d &torque)
{
    int id1 = obj1->n6DOF; int id2 = obj2->n6DOF;
    Eigen::Vector3d center1 = obj1->c_; Eigen::Vector3d center2 = obj2->c_;
    Eigen::Vector3d r1 = contact_point - center1; Eigen::Vector3d r2 = contact_point - center2;
    Eigen::Vector3d v1(obj1->p_(0)/obj1->Mass_fb, obj1->p_(1)/obj1->Mass_fb, obj1->p_(2)/obj1->Mass_fb);
    Eigen::Vector3d v2(obj2->p_(0)/obj2->Mass_fb, obj2->p_(1)/obj2->Mass_fb, obj2->p_(2)/obj2->Mass_fb);
    Eigen::Vector3d omega1 = obj1->omega_I; Eigen::Vector3d omega2 = obj2->omega_I;
    Eigen::Vector3d v1_contact = v1 + omega1.cross(r1); Eigen::Vector3d v2_contact = v2 + omega2.cross(r2);
    Eigen::Vector3d v_rel = v2_contact - v1_contact; double v_rel_n = v_rel.dot(normal);
    Eigen::Vector3d v_rel_t = v_rel - v_rel_n * normal;

    auto key = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    if(contact_history.find(key) == contact_history.end()) { ContactHistory h; h.tangential_overlap.setZero(); h.in_contact=true; h.last_update_time=p->simtime; contact_history[key]=h; }
    double dt = p->simtime - contact_history[key].last_update_time; if(dt<=0.0) dt=p->dt;
    contact_history[key].in_contact=true; contact_history[key].last_update_time=p->simtime;

    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);
    double G = young_modulus / (1.0 + poisson_ratio);
    //double G = young_modulus / (2.0 * (1.0 + poisson_ratio));
    double mi = obj1->Mass_fb; double mj = obj2->Mass_fb; double meff = (mi*mj)/(mi+mj);
    double K_hertz = (4.0/3.0) * E_eff * std::sqrt(R_eff);
    const double en = std::min(std::max(restitution_coefficient, 1.0e-6), 0.999999);
    double ethan = -2.2664 * std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));
    Eigen::Vector3d Fn = (-(4.0/3.0) * E_eff * std::sqrt(R_eff) * std::pow(overlap, 1.5)
                          - std::sqrt(meff * K_hertz) * ethan * std::pow(overlap, 0.25) * v_rel_n) * normal;

    auto &overlap_t = contact_history[key].tangential_overlap;
    if(v_rel_t.norm() > 1.0e-10)
    {
        overlap_t += v_rel_t * dt; overlap_t -= normal * normal.dot(overlap_t);
        double kt = 8.0 * G * std::sqrt(R_eff * overlap);
        Eigen::Vector3d Ft = - kt * overlap_t; // no damping
        double ft = Ft.norm(); double ft_max = friction_coefficient * Fn.norm();
        if(ft > ft_max)
        {
            if(overlap_t.norm() > 0.0) { Ft *= (ft_max/ft); }
            else { Ft.setZero(); }
        }
        force = Fn + Ft;
    }
    else { force = Fn; }
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
    
    // Calculate effective radius and effective Young's modulus (identical materials assumed)
    const double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    const double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);

    // Hertz contact elastic term
    const double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    const double fn_elastic = (4.0/3.0) * k_hertz * pow(overlap, 1.5);

    // Mindlin-style normal damping consistent with LIGGGHTS
    // Sn = 2*E_eff*sqrt(R_eff*delta)
    const double Sn = 2.0 * E_eff * sqrt(R_eff * overlap);
    const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
    // beta from restitution coefficient (clamped)
    const double en = std::min(std::max(restitution_coefficient, 1.0e-6), 0.999999);
    const double beta = std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));
    const double sqrtFiveOverSix = 0.9128709291752769;
    const double gamman = -2.0 * sqrtFiveOverSix * beta * std::sqrt(Sn * meff);
    const double fn_damping = -gamman * v_rel_n;
    double fn = fn_elastic + fn_damping;
    fn = std::max(fn, 0.0); // Ensure normal force is repulsive
    
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
        // Increment and project tangential overlap into tangential plane
        auto &th = contact_history[std::make_pair(std::min(id1, id2), std::max(id1, id2))].tangential_overlap;
        th += v_rel_t * dt;
        th -= normal * normal.dot(th);

        // Mindlin tangential stiffness and damping
        // G = E/(2(1+nu)); for identical materials: G* = G / (2(2 - nu))
        const double G = young_modulus / (2.0 * (1.0 + poisson_ratio));
        const double Geff_mindlin = G / (2.0 * (2.0 - poisson_ratio));
        const double kt = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
        const double St = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
        const double gammat = -2.0 * sqrtFiveOverSix * beta * std::sqrt(St * meff);

        Eigen::Vector3d ft_vector = kt * th - gammat * v_rel_t;

        // Coulomb limit
        const double ft_mag = ft_vector.norm();
        const double ft_max = friction_coefficient * fn;
        if(ft_mag > ft_max)
        {
            ft_vector *= (ft_max / ft_mag);
            th *= (ft_max / ft_mag);
        }

        // Total force
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
    
    // Calculate effective radius and moduli (identical materials assumed)
    const double R_eff = calculate_effective_radius(obj1->radius, obj2->radius);
    const double E_eff = calculate_effective_young_modulus(young_modulus, young_modulus, poisson_ratio, poisson_ratio);

    // Hertz normal component with LIGGGHTS-like damping
    const double k_hertz = calculate_hertz_stiffness(E_eff, R_eff);
    const double fn_elastic = (4.0/3.0) * k_hertz * pow(overlap, 1.5);
    const double Sn = 2.0 * E_eff * std::sqrt(R_eff * overlap);
    const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
    const double en = std::min(std::max(restitution_coefficient, 1.0e-6), 0.999999);
    const double beta = std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));
    const double sqrtFiveOverSix = 0.9128709291752769;
    const double gamman = -2.0 * sqrtFiveOverSix * beta * std::sqrt(Sn * meff);
    const double fn_damping = -gamman * v_rel_n;
    double fn = fn_elastic + fn_damping;
    fn = std::max(fn, 0.0);

    // Mindlin tangential stiffness and damping
    const double G = young_modulus / (2.0 * (1.0 + poisson_ratio));
    const double Geff_mindlin = G / (2.0 * (2.0 - poisson_ratio));
    const double k_t = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
    const double St = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
    const double gammat = -2.0 * sqrtFiveOverSix * beta * std::sqrt(St * meff);
    
    // Update tangential overlap (displacement)
    if(v_rel_t_mag > 1.0e-10)
    {
        Eigen::Vector3d t_hat = v_rel_t / v_rel_t_mag;
        
        // Increment tangential overlap based on relative velocity
        contact_history[contact_pair].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[contact_pair].tangential_overlap -= 
            normal * normal.dot(contact_history[contact_pair].tangential_overlap);
        
        // Calculate the tangential force based on tangential spring plus viscous term
        Eigen::Vector3d ft_vector = k_t * contact_history[contact_pair].tangential_overlap - gammat * v_rel_t;
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
    // Calculate base Hertzian stiffness factor (without 4/3): k = E_eff * sqrt(R_eff)
    // The (4/3) factor is applied at the point of computing the normal force.
    return E_eff * sqrt(R_eff);
}

// NEW: PacIFiC-style effective properties calculation
double sixdof_collision::calculate_effective_shear_modulus(double E1, double E2, double nu1, double nu2)
{
    // Calculate effective shear modulus based on PacIFiC implementation
    // G* = 1 / (2(2-ν₀)(1+ν₀)/E₀ + 2(2-ν₁)(1+ν₁)/E₁)
    double G1 = E1 / (2.0 * (1.0 + nu1));
    double G2 = E2 / (2.0 * (1.0 + nu2));
    
    // For most materials, G* ≈ 4*E* (as noted in PacIFiC documentation)
    double E_eff = calculate_effective_young_modulus(E1, E2, nu1, nu2);
    return 4.0 * E_eff;
}

double sixdof_collision::calculate_effective_mass(double m1, double m2)
{
    // Calculate effective mass: m* = 1 / (1/m₁ + 1/m₂)
    return 1.0 / (1.0/m1 + 1.0/m2);
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
    
    // BVH-accelerated path: if both objects have BVH, use it for fast rejection
    if(obj1->use_bvh && obj1->mesh_bvh && obj1->mesh_bvh->is_built())
    {
        // Check if obj2's sphere intersects obj1's BVH
        if(!obj1->mesh_bvh->intersects_sphere(obj2->c_, obj2->radius))
        {
            return false;  // BVH says no collision possible
        }
    }
    
    if(obj2->use_bvh && obj2->mesh_bvh && obj2->mesh_bvh->is_built())
    {
        // Check if obj1's sphere intersects obj2's BVH
        if(!obj2->mesh_bvh->intersects_sphere(obj1->c_, obj1->radius))
        {
            return false;  // BVH says no collision possible
        }
    }
        
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
    // INDUSTRY-STANDARD TRIANGLE-TRIANGLE INTERSECTION
    // Based on: Ericson "Real-Time Collision Detection" (2005), Chapter 5.3.4
    // Uses Separating Axis Theorem (SAT) with 11 potential separating axes
    //
    // Theory: Two convex objects don't intersect IFF there exists a separating axis
    //         where their projections don't overlap.
    //
    // For triangles, potential separating axes are:
    //   1-2:  Face normals of both triangles (2 axes)
    //   3-11: Cross products of all edge pairs (3×3 = 9 axes)
    //
    // This is the same algorithm used in: PhysX, Bullet, Unity, Unreal Engine
    
    const double EPSILON = 1e-8;
    
    // Calculate triangle normals
    Eigen::Vector3d n1 = (v1[1] - v1[0]).cross(v1[2] - v1[0]);
    Eigen::Vector3d n2 = (v2[1] - v2[0]).cross(v2[2] - v2[0]);
    
    // Check for degenerate triangles (zero area)
    double n1_len = n1.norm();
    double n2_len = n2.norm();
    if(n1_len < EPSILON || n2_len < EPSILON)
        return false;  // Degenerate triangle
    
    n1 /= n1_len;  // Normalize
    n2 /= n2_len;
    
    // Calculate triangle edges
    Eigen::Vector3d e1[3], e2[3];
    e1[0] = v1[1] - v1[0];
    e1[1] = v1[2] - v1[1];
    e1[2] = v1[0] - v1[2];
    
    e2[0] = v2[1] - v2[0];
    e2[1] = v2[2] - v2[1];
    e2[2] = v2[0] - v2[2];
    
    // Track minimum penetration depth and corresponding axis
    double min_penetration = 1e10;
    Eigen::Vector3d best_axis;
    bool found_separation = false;
    
    // Helper lambda: Test separation along an axis
    auto test_axis = [&](const Eigen::Vector3d& axis) -> bool
    {
        double axis_len = axis.norm();
        if(axis_len < EPSILON)
            return true;  // Degenerate axis, skip
        
        Eigen::Vector3d axis_norm = axis / axis_len;
        
        // Project triangle 1 vertices onto axis
        double proj1[3];
        for(int i = 0; i < 3; ++i)
            proj1[i] = v1[i].dot(axis_norm);
        
        double min1 = std::min({proj1[0], proj1[1], proj1[2]});
        double max1 = std::max({proj1[0], proj1[1], proj1[2]});
        
        // Project triangle 2 vertices onto axis
        double proj2[3];
        for(int i = 0; i < 3; ++i)
            proj2[i] = v2[i].dot(axis_norm);
        
        double min2 = std::min({proj2[0], proj2[1], proj2[2]});
        double max2 = std::max({proj2[0], proj2[1], proj2[2]});
        
        // Check for separation
        if(max1 < min2 - EPSILON || max2 < min1 - EPSILON)
        {
            found_separation = true;
            return false;  // Separating axis found - no collision
        }
        
        // Calculate penetration depth along this axis
        double penetration = std::min(max1 - min2, max2 - min1);
        
        // Track axis with minimum penetration (this gives us the best contact normal)
        if(penetration < min_penetration)
        {
            min_penetration = penetration;
            best_axis = axis_norm;
            
            // Make sure normal points from triangle 2 to triangle 1
            Eigen::Vector3d center1 = (v1[0] + v1[1] + v1[2]) / 3.0;
            Eigen::Vector3d center2 = (v2[0] + v2[1] + v2[2]) / 3.0;
            if((center1 - center2).dot(best_axis) < 0)
                best_axis = -best_axis;
        }
        
        return true;  // Continue testing
    };
    
    // TEST 1-2: Face normals (2 axes)
    if(!test_axis(n1)) return false;
    if(!test_axis(n2)) return false;
    
    // TEST 3-11: Edge-edge cross products (9 axes)
    // For each edge of triangle 1 crossed with each edge of triangle 2
    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            Eigen::Vector3d axis = e1[i].cross(e2[j]);
            if(!test_axis(axis)) return false;
        }
    }
    
    // If we get here, no separating axis was found → triangles intersect!
    
    // Set output parameters
    overlap = min_penetration;
    normal = best_axis;
    
    // Calculate contact point as the average of the closest points
    // For simplicity, use the centroid of the overlapping region
    // (More sophisticated methods exist, but this is sufficient for collision response)
    
    // Find vertices of triangle 1 that are "inside" triangle 2 (based on normal)
    std::vector<Eigen::Vector3d> contact_points;
    
    for(int i = 0; i < 3; ++i)
    {
        // Simple heuristic: check if vertex is close to the other triangle's plane
        double dist = (v1[i] - v2[0]).dot(n2);
        if(fabs(dist) < min_penetration + EPSILON)
            contact_points.push_back(v1[i]);
        
        dist = (v2[i] - v1[0]).dot(n1);
        if(fabs(dist) < min_penetration + EPSILON)
            contact_points.push_back(v2[i]);
    }
    
    // Calculate contact point
    if(contact_points.empty())
    {
        // Fallback: use midpoint of triangle centroids
        Eigen::Vector3d center1 = (v1[0] + v1[1] + v1[2]) / 3.0;
        Eigen::Vector3d center2 = (v2[0] + v2[1] + v2[2]) / 3.0;
        contact = 0.5 * (center1 + center2);
    }
    else
    {
        // Average of contact points
        contact = Eigen::Vector3d::Zero();
        for(const auto& pt : contact_points)
            contact += pt;
        contact /= contact_points.size();
    }
    
    return true;  // Collision detected!
}

bool sixdof_collision::detect_collision_adaptive(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                 Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // HYBRID ADAPTIVE COLLISION DETECTION
    // Automatically selects the best algorithm based on object complexity
    
    // Get triangle counts for both objects
    int tri_count_1 = obj1->tricount;
    int tri_count_2 = obj2->tricount;
    int max_tri_count = std::max(tri_count_1, tri_count_2);
    
    // Strategy selection based on complexity
    if(max_tri_count < collision_simple_threshold)
    {
        // CASE 1: Both objects are simple (< 50 triangles)
        // Use fast sphere-sphere approximation
        // This is very fast and good enough for simple geometries
        
        return detect_collision(p, pgc, obj1, obj2, contact_point, normal, overlap);
    }
    else if(max_tri_count < collision_moderate_threshold)
    {
        // CASE 2: Moderate complexity (50-500 triangles)
        // Use triangle-triangle with BVH acceleration
        // BVH makes this fast enough for real-time
        
        return detect_triangle_collision(p, pgc, obj1, obj2, contact_point, normal, overlap);
    }
    else
    {
        // CASE 3: High complexity (> 500 triangles)
        // Use triangle-triangle with BVH acceleration
        // BVH is essential here - without it this would be too slow
        
        return detect_triangle_collision(p, pgc, obj1, obj2, contact_point, normal, overlap);
    }
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

// NEW: Enhanced PacIFiC-based Hertz contact force model
void sixdof_collision::calculate_pacific_hertz_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // PacIFiC-style effective properties calculation
    double Req = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double avmass = calculate_effective_mass(obj1->Mass_fb, obj2->Mass_fb);
    double deltan = overlap;  // Overlap distance (positive for contact)
    double sqrtReqdeltan = sqrt(Req * deltan);
    
    // PacIFiC-style stiffness parameters
    double Sn = 2.0 * pacific_Es * sqrtReqdeltan;  // Normal stiffness parameter
    double St = 8.0 * pacific_Gs * sqrtReqdeltan;  // Tangential stiffness parameter
    
    // Normal non-linear elastic force (PacIFiC style)
    double kn = (4.0 / 3.0) * pacific_Es * sqrtReqdeltan;
    Eigen::Vector3d delFN = kn * overlap * normal;
    
    // Normal non-linear dissipative force (PacIFiC style)
    double gamman = pacific_m2sqrt56 * pacific_beta * sqrt(avmass * Sn);
    delFN -= gamman * v_rel_n * normal;
    
    double normFN = delFN.norm();
    
    // Tangential non-linear dissipative force (PacIFiC style)
    double gammat = pacific_m2sqrt56 * pacific_beta * sqrt(avmass * St);
    Eigen::Vector3d delFT = -gammat * v_rel_t;
    
    // Tangential Coulomb saturation (PacIFiC style)
    double fn = pacific_muc * normFN;
    double ft = delFT.norm();
    
    if (fn < ft && v_rel_t_mag > 1.0e-10) {
        // Unit tangential vector in the reverse direction of the relative velocity
        Eigen::Vector3d tangent = -v_rel_t / v_rel_t_mag;
        delFT = tangent * fn;
    }
    
    // Total force vector
    force = delFN + delFT;
    
    // PacIFiC-style rolling resistance torque
    if (pacific_kr > 0.0) {
        // Relative angular velocity at contact point
        Eigen::Vector3d wrel = omega1 - omega2;
        double normwrel = wrel.norm();
        
        // Tangential velocity contribution from angular velocities
        Eigen::Vector3d wt1 = omega1.cross(r1);
        Eigen::Vector3d wt2 = omega2.cross(r2);
        Eigen::Vector3d wtrel = wt1 - wt2;
        double normwtrel = wtrel.norm();
        
        if (normwrel > 1.0e-10) {
            Eigen::Vector3d delM = -(pacific_kr * Req * normFN * normwtrel / normwrel) * wrel;
            torque = r1.cross(force) + delM;
        } else {
            torque = r1.cross(force);
        }
    } else {
        torque = r1.cross(force);
    }
    
    // Update contact history for tangential forces
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end()) {
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        history.previous_normal = normal;
        history.tangential_spring.setZero();
        history.contact_duration = 0.0;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    } else {
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].contact_duration += p->dt;
    }
}

// NEW: Enhanced PacIFiC-based Hooke contact force model
void sixdof_collision::calculate_pacific_hooke_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // PacIFiC-style effective properties calculation
    double avmass = calculate_effective_mass(obj1->Mass_fb, obj2->Mass_fb);
    
    // Normal linear elastic force (PacIFiC Hooke style)
    Eigen::Vector3d delFN = hooke_kn * overlap * normal;
    
    // Normal dissipative force (derive damping ratio from restitution)
    // zeta = -ln(e) / sqrt(pi^2 + ln^2(e)); c = 2*zeta*sqrt(k*meff)
    const double en = std::min(std::max(hooke_en, 1.0e-6), 0.999999);
    const double zeta = -std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));
    double gamman = 2.0 * zeta * std::sqrt(hooke_kn * avmass);
    delFN -= gamman * v_rel_n * normal;
    
    double normFN = delFN.norm();
    
    // Get or create contact history for this pair
    auto it = contact_history.find(std::make_pair(min(id1, id2), max(id1, id2)));
    if(it == contact_history.end()) {
        ContactHistory history;
        history.tangential_overlap.setZero();
        history.in_contact = true;
        history.last_update_time = p->simtime;
        history.previous_normal = normal;
        history.tangential_spring.setZero();
        history.contact_duration = 0.0;
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))] = history;
    }
    
    // Get time step
    double dt = p->simtime - contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time;
    if(dt <= 0.0) dt = p->dt;
    
    // Update contact history
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].in_contact = true;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].last_update_time = p->simtime;
    contact_history[std::make_pair(min(id1, id2), max(id1, id2))].contact_duration += p->dt;
    
    // Tangential force calculation with memory effects (PacIFiC Hooke style)
    Eigen::Vector3d delFT;
    
    if(v_rel_t_mag > 1.0e-10) {
        // Update tangential overlap based on relative velocity
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap += v_rel_t * dt;
        
        // Project tangential overlap to the current tangential plane
        contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap -= 
            normal * normal.dot(contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap);
        
        // Calculate the tangential force based on tangential spring
        Eigen::Vector3d ft_spring = hooke_kt * contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap;
        
        // Calculate the tangential damping force
        // Tangential damping from a tangential restitution or by mirroring normal zeta
        double zeta_t;
        if (hooke_et == -1.0) {
            zeta_t = zeta; // use same damping ratio if not provided
        } else {
            const double et = std::min(std::max(hooke_et, 1.0e-6), 0.999999);
            zeta_t = -std::log(et) / std::sqrt(M_PI*M_PI + std::log(et)*std::log(et));
        }
        double gammat = 2.0 * zeta_t * std::sqrt(hooke_kt * avmass);
        Eigen::Vector3d ft_damp = -gammat * v_rel_t;
        
        delFT = ft_spring + ft_damp;
        
        // Tangential Coulomb saturation
        double ft_mag = delFT.norm();
        double ft_max = hooke_muc * normFN;
        
        if(ft_mag > ft_max) {
            // Scale tangential force to the maximum allowed
            delFT *= (ft_max / ft_mag);
            
            // Update tangential overlap to match the maximum force
            contact_history[std::make_pair(min(id1, id2), max(id1, id2))].tangential_overlap *= (ft_max / ft_mag);
        }
    } else {
        delFT.setZero();
    }
    
    // Total force vector
    force = delFN + delFT;
    
    // PacIFiC-style rolling resistance torque
    if (hooke_kr > 0.0) {
        // Relative angular velocity at contact point
        Eigen::Vector3d wrel = omega1 - omega2;
        double normwrel = wrel.norm();
        
        // Tangential velocity contribution from angular velocities
        Eigen::Vector3d wt1 = omega1.cross(r1);
        Eigen::Vector3d wt2 = omega2.cross(r2);
        Eigen::Vector3d wtrel = wt1 - wt2;
        double normwtrel = wtrel.norm();
        
        if (normwrel > 1.0e-10) {
            double Req = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
            Eigen::Vector3d delM = -(hooke_kr * Req * normFN * normwtrel / normwrel) * wrel;
            torque = r1.cross(force) + delM;
        } else {
            torque = r1.cross(force);
        }
    } else {
        torque = r1.cross(force);
    }
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

// NEW: Clear collision forces and torques for all objects
void sixdof_collision::clear_collision_forces()
{
    for(int i = 0; i < max_objects; ++i)
    {
        collision_forces[i].setZero();
        collision_torques[i].setZero();
    }
}

// NEW: Broadcast collision forces and torques from rank 0 to all processors
void sixdof_collision::broadcast_collision_forces(lexer *p, ghostcell *pgc)
{
    // Broadcast collision forces and torques for all objects
    for(int i = 0; i < max_objects; ++i)
    {
        // Broadcast forces (3 components: x, y, z)
        MPI_Bcast(&collision_forces[i](0), 3, MPI_DOUBLE, 0, pgc->mpi_comm);
        
        // Broadcast torques (3 components: x, y, z)
        MPI_Bcast(&collision_torques[i](0), 3, MPI_DOUBLE, 0, pgc->mpi_comm);
    }
    
    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"6DOF Collision: Broadcasted forces and torques to all processors"<<endl;
    }
}

// NEW: Verify that all processors have synchronized collision forces (debug function)
void sixdof_collision::verify_collision_forces_synchronization(lexer *p, ghostcell *pgc)
{
    // This function can be called to verify that MPI communication worked correctly
    // It's useful for debugging but not needed for normal operation
    
    bool forces_synchronized = true;
    
    for(int i = 0; i < max_objects; ++i)
    {
        // Check if forces are synchronized across all processors
        double force_x_sum = 0.0, force_y_sum = 0.0, force_z_sum = 0.0;
        double torque_x_sum = 0.0, torque_y_sum = 0.0, torque_z_sum = 0.0;
        
        // Sum forces across all processors
        MPI_Allreduce(&collision_forces[i](0), &force_x_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&collision_forces[i](1), &force_y_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&collision_forces[i](2), &force_z_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        
        MPI_Allreduce(&collision_torques[i](0), &torque_x_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&collision_torques[i](1), &torque_y_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&collision_torques[i](2), &torque_z_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        
        // Check if the sum equals the value times number of processors
        int num_procs;
        MPI_Comm_size(pgc->mpi_comm, &num_procs);
        
        double expected_force_x = collision_forces[i](0) * num_procs;
        double expected_force_y = collision_forces[i](1) * num_procs;
        double expected_force_z = collision_forces[i](2) * num_procs;
        
        double expected_torque_x = collision_torques[i](0) * num_procs;
        double expected_torque_y = collision_torques[i](1) * num_procs;
        double expected_torque_z = collision_torques[i](2) * num_procs;
        
        // Check for synchronization (allow small numerical tolerance)
        double tolerance = 1.0e-10;
        if(fabs(force_x_sum - expected_force_x) > tolerance || 
           fabs(force_y_sum - expected_force_y) > tolerance || 
           fabs(force_z_sum - expected_force_z) > tolerance ||
           fabs(torque_x_sum - expected_torque_x) > tolerance || 
           fabs(torque_y_sum - expected_torque_y) > tolerance || 
           fabs(torque_z_sum - expected_torque_z) > tolerance)
        {
            forces_synchronized = false;
        }
    }
    
    if(p->mpirank == 0)
    {
        if(forces_synchronized)
        {
            cout<<"6DOF Collision: SUCCESS - All collision forces are synchronized across processors"<<endl;
        }
        else
        {
            cout<<"6DOF Collision: ERROR - Collision forces are NOT synchronized across processors!"<<endl;
        }
    }
} 