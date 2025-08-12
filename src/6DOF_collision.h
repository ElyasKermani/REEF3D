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
class sixdof_collision_grid;

using namespace std;

// Enumeration for different collision models
enum class ContactForceModel {
    Linear,      // Linear spring-dashpot model
    Hertz,       // Non-linear Hertzian elastic contact
    HertzMindlin, // Hertz with tangential history
    DMT,         // Derjaguin-Muller-Toporov model for adhesive contacts
    JKR,         // Johnson-Kendall-Roberts model for strong adhesion
    PacIFiCHertz, // Enhanced Hertz model based on PacIFiC implementation
    PacIFiCHooke  // Enhanced Hooke model based on PacIFiC implementation
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
    
private:

    // Detect collision between two 6DOF objects
    bool detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2, 
                         Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);
    
    // Detect collision using triangle mesh data
    bool detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                 Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);
    
    // Check for intersection between two triangles
    bool triangle_triangle_intersection(const Eigen::Vector3d v1[3], const Eigen::Vector3d v2[3],
                                      Eigen::Vector3d &contact, Eigen::Vector3d &normal, double &overlap);
    
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
    
    // Calculate JKR contact force (strong adhesion)
    void calculate_jkr_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                   const Eigen::Vector3d &contact_point, 
                                   const Eigen::Vector3d &normal, 
                                   const double overlap,
                                   Eigen::Vector3d &force, 
                                   Eigen::Vector3d &torque);
    
    // NEW: Enhanced PacIFiC-based Hertz contact force model
    void calculate_pacific_hertz_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                             const Eigen::Vector3d &contact_point, 
                                             const Eigen::Vector3d &normal, 
                                             const double overlap,
                                             Eigen::Vector3d &force, 
                                             Eigen::Vector3d &torque);
    
    // NEW: Enhanced PacIFiC-based Hooke contact force model
    void calculate_pacific_hooke_contact_force(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // NEW: PacIFiC-style effective properties calculation
    double calculate_effective_shear_modulus(double E1, double E2, double nu1, double nu2);
    double calculate_effective_mass(double m1, double m2);
    
    // Contact model parameters
    ContactForceModel contact_model;
    
    // Linear model parameters
    double spring_constant_n;        // Normal spring constant
    double spring_constant_t;        // Tangential spring constant
    double damping_constant_n;       // Normal damping constant
    double damping_constant_t;       // Tangential damping constant
    double friction_coefficient;     // Friction coefficient
    double restitution_coefficient;  // Restitution coefficient
    
    // Rolling friction parameters
    double rolling_friction_coefficient;  // Rolling friction coefficient
    double rolling_stiffness;            // Rolling spring constant
    double rolling_damping;              // Rolling damping constant
    double rolling_torque_threshold;     // Threshold for rolling torque activation
    
    // Hertz model parameters
    double young_modulus;            // Young's modulus
    double poisson_ratio;            // Poisson's ratio
    
    // DMT model parameters
    double surface_energy;           // Surface energy
    double dmt_cutoff_threshold;     // Cutoff threshold for DMT
    
    // JKR model parameters
    double surface_energy_jkr;     // Surface energy for JKR model
    double jkr_cutoff_threshold;   // Cutoff threshold for JKR model
    
    // NEW: PacIFiC-style enhanced parameters
    double pacific_Es;              // Effective Young's modulus (PacIFiC style)
    double pacific_en;              // Normal restitution coefficient
    double pacific_Gs;              // Effective shear modulus
    double pacific_muc;             // Coulomb friction coefficient
    double pacific_kr;              // Rolling resistance coefficient
    double pacific_beta;            // Damping factor: log(en)/sqrt(π² + log²(en))
    double pacific_m2sqrt56;        // Constant: -2*sqrt(5/6)
    
    // Hooke model parameters (PacIFiC style)
    double hooke_kn;                // Normal stiffness coefficient
    double hooke_kt;                // Tangential stiffness coefficient
    double hooke_en;                // Normal restitution coefficient
    double hooke_et;                // Tangential restitution coefficient
    double hooke_muc;               // Coulomb friction coefficient
    double hooke_kr;                // Rolling resistance coefficient
    
    // Sub-stepping parameters
    bool use_substeps;
    int max_substeps;
    
    // NEW: Collision force synchronization method
    enum class SyncMethod {
        GlobalSum,      // Use pgc->globalsum (most robust, slower)
        Broadcast       // Use MPI_Bcast from rank 0 (faster, requires rank 0 calculation)
    };
    SyncMethod sync_method;
    
    // Contact history for tangential forces
    struct ContactHistory {
        Eigen::Vector3d tangential_overlap;
        bool in_contact;
        double last_update_time;
        // NEW: Enhanced contact history for PacIFiC-style models
        Eigen::Vector3d previous_normal;      // Previous contact normal
        Eigen::Vector3d tangential_spring;    // Tangential spring displacement
        double contact_duration;              // Duration of contact
    };
    map<pair<int, int>, ContactHistory> contact_history;
    
    // Clear contact history for pairs no longer in contact
    void update_contact_history(lexer *p);
    
    // Grid-based collision detection system
    sixdof_collision_grid* collision_grid;
    
    // For distance calculation
    double calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2);
    
    // Sub-time stepping for resolving collisions at smaller timesteps
    void resolve_collision_with_substeps(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                       const Eigen::Vector3d &contact_point, 
                                       const Eigen::Vector3d &normal, 
                                       const double overlap,
                                       Eigen::Vector3d &force, 
                                       Eigen::Vector3d &torque);
    
    // Velocity-Verlet integration step for collision resolution
    void velocity_verlet_step(lexer *p, ghostcell *pgc, sixdof_obj *obj, 
                            const Eigen::Vector3d &force, 
                            const Eigen::Vector3d &torque, 
                            double dt);
    
    // Calculate rolling friction torque
    void calculate_rolling_friction_torque(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                         const Eigen::Vector3d &contact_point,
                                         const Eigen::Vector3d &normal,
                                         const double overlap,
                                         Eigen::Vector3d &rolling_torque);
                                         
    // Calculate twisting resistance
    void calculate_twisting_resistance(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                     const Eigen::Vector3d &contact_point,
                                     const Eigen::Vector3d &normal,
                                     const double overlap,
                                     Eigen::Vector3d &twisting_torque);
    
    // NEW: PacIFiC-style rolling resistance calculation
    void calculate_pacific_rolling_resistance(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                            const Eigen::Vector3d &contact_point,
                                            const Eigen::Vector3d &normal,
                                            const double overlap,
                                            const double normFN,
                                            Eigen::Vector3d &rolling_torque);
    
    // NEW: Function to synchronize collision forces across all processors
    void synchronize_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
    // NEW: Set synchronization method for performance tuning
    void set_sync_method(SyncMethod method) { sync_method = method; }
    
    // NEW: Get current synchronization method
    SyncMethod get_sync_method() const { return sync_method; }
    
    // Set contact force model
    void set_contact_force_model(ContactForceModel model) { contact_model = model; }
    
    // Get current contact force model
    ContactForceModel get_contact_force_model() const { return contact_model; }
};

#endif 