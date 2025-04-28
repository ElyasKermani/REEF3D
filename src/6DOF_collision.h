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
#include<map>

class lexer;
class ghostcell;
class sixdof_obj;

using namespace std;

// Enumeration for different collision models
enum class ContactForceModel {
    Linear,      // Linear spring-dashpot model
    Hertz,       // Non-linear Hertzian elastic contact
    HertzMindlin, // Hertz with tangential history
    DMT         // Derjaguin-Muller-Toporov model for adhesive contacts
};

class sixdof_collision
{
public:

    sixdof_collision(lexer *p, ghostcell *pgc);
    virtual ~sixdof_collision();
    
    // Calculate collision forces between all 6DOF objects
    void calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
    // Set the contact force model to use
    void set_contact_force_model(ContactForceModel model) { contact_model = model; }
    
    // New function for size-dependent parameter calculation
    void calculate_size_dependent_parameters(lexer *p, sixdof_obj *obj1, sixdof_obj *obj2);
    
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
    
    // Calculate Hertzian contact force (non-linear elastic)
    void calculate_hertz_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                     const Eigen::Vector3d &contact_point, 
                                     const Eigen::Vector3d &normal, 
                                     const double overlap,
                                     Eigen::Vector3d &force, 
                                     Eigen::Vector3d &torque);
    
    // Calculate Hertz-Mindlin contact force with tangential history
    void calculate_hertz_mindlin_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                            const Eigen::Vector3d &contact_point, 
                                            const Eigen::Vector3d &normal, 
                                            const double overlap,
                                            Eigen::Vector3d &force, 
                                            Eigen::Vector3d &torque);
    
    // Calculate DMT (Derjaguin-Muller-Toporov) model for adhesive contacts
    void calculate_dmt_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                   const Eigen::Vector3d &contact_point, 
                                   const Eigen::Vector3d &normal, 
                                   const double overlap,
                                   Eigen::Vector3d &force, 
                                   Eigen::Vector3d &torque);
    
    // Calculate effective material properties
    double calculate_effective_young_modulus(double E1, double E2, double nu1, double nu2);
    double calculate_effective_radius(double R1, double R2);
    
    // Helper function for Hertzian contact
    double calculate_hertz_stiffness(double E_eff, double R_eff);
    
    // Parameters for the collision models
    ContactForceModel contact_model;
    
    // Base parameters for scaling
    double base_spring_constant;
    double base_damping_constant;
    double base_young_modulus;
    double base_surface_energy;
    
    // Current parameters (scaled)
    double scaled_spring_constant;
    double scaled_damping_constant;
    double scaled_young_modulus;
    double scaled_surface_energy;
    double effective_young_modulus;
    double effective_radius;
    
    // Fixed parameters
    double poisson_ratio;
    double friction_coefficient;
    double restitution_coefficient;
    double dmt_cutoff_threshold;
    
    // Rolling resistance parameters
    double rolling_friction_coefficient;
    double rolling_viscous_damping;
    double scaled_rolling_damping;
    
    // Contact history for Hertz-Mindlin model
    struct ContactHistory {
        Eigen::Vector3d tangential_overlap;
        bool in_contact;
        double last_update_time;
    };
    
    // Map to store contact history between object pairs
    map<pair<int, int>, ContactHistory> contact_history;
    
    // Clear contact history for pairs no longer in contact
    void update_contact_history(lexer *p);
    
    // Simplified bounding sphere collision detection
    double calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2);
    
    // For simplified spherical collision (using bounding spheres)
    vector<double> bounding_radius;
};

#endif 