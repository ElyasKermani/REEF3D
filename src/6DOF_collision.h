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

// AABB (Axis-Aligned Bounding Box) structure for fast collision pre-filtering
struct AABB {
    Eigen::Vector3d min;  // Minimum corner (x_min, y_min, z_min)
    Eigen::Vector3d max;  // Maximum corner (x_max, y_max, z_max)
    
    // Default constructor
    AABB() : min(0,0,0), max(0,0,0) {}
    
    // Update AABB from center and radius (sphere)
    inline void update_from_sphere(const Eigen::Vector3d& center, double radius) {
        min = center - Eigen::Vector3d(radius, radius, radius);
        max = center + Eigen::Vector3d(radius, radius, radius);
    }
    
    // Check if two AABBs overlap (fast early rejection test)
    inline bool overlaps(const AABB& other) const {
        return (max.x() >= other.min.x() && min.x() <= other.max.x()) &&
               (max.y() >= other.min.y() && min.y() <= other.max.y()) &&
               (max.z() >= other.min.z() && min.z() <= other.max.z());
    }
    
    // Expand AABB by a margin (useful for tolerance)
    inline void expand(double margin) {
        min.x() -= margin;
        min.y() -= margin;
        min.z() -= margin;
        max.x() += margin;
        max.y() += margin;
        max.z() += margin;
    }
};

// Enumeration for different collision models
enum class ContactForceModel {
    Linear,      // Linear spring-dashpot model
    Hertz,       // Non-linear Hertzian elastic contact
    HertzMindlin, // Hertz with tangential history
    PhasicFlowLinearLimited,      // phasicFlow linear model with history limiting
    PhasicFlowLinearNonLimited,   // phasicFlow linear model without history limiting
    PhasicFlowNonLinearLimited,   // phasicFlow Hertz–Mindlin with history limiting (no tangential damping)
    PhasicFlowNonLinearNonLimited // phasicFlow Hertz–Mindlin without history limiting
};

class sixdof_collision
{
public:

    sixdof_collision(lexer *p, ghostcell *pgc);
    virtual ~sixdof_collision();
    
    // Calculate collision forces between all 6DOF objects
    void calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
    // NEW: Calculate ground contact forces for all 6DOF objects
    void calculate_ground_contact_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
    // NEW: Calculate boundary wall contact forces (side walls)
    void calculate_boundary_wall_contact_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);
    
    // Set the contact force model to use
    void set_contact_force_model(ContactForceModel model) { contact_model = model; }
    
private:

    // Detect collision between two 6DOF objects
    bool detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2, 
                         Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);
    
    // Detect collision using triangle mesh data
    bool detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                 Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);
    
    // Adaptive collision detection - chooses algorithm based on object complexity
    bool detect_collision_adaptive(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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

    // phasicFlow-equivalent models
    void calculate_phasicflow_linear_limited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                           const Eigen::Vector3d &contact_point,
                                           const Eigen::Vector3d &normal,
                                           const double overlap,
                                           Eigen::Vector3d &force,
                                           Eigen::Vector3d &torque);

    void calculate_phasicflow_linear_nonlimited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                              const Eigen::Vector3d &contact_point,
                                              const Eigen::Vector3d &normal,
                                              const double overlap,
                                              Eigen::Vector3d &force,
                                              Eigen::Vector3d &torque);

    void calculate_phasicflow_nonlinear_limited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                              const Eigen::Vector3d &contact_point,
                                              const Eigen::Vector3d &normal,
                                              const double overlap,
                                              Eigen::Vector3d &force,
                                              Eigen::Vector3d &torque);

    void calculate_phasicflow_nonlinear_nonlimited(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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
    
    // Contact model parameters
    ContactForceModel contact_model;
    
    // Linear model parameters
    double spring_constant_n;        // Normal spring constant
    double spring_constant_t;        // Tangential spring constant
    double damping_constant_n;       // Normal damping constant
    double damping_constant_t;       // Tangential damping constant
    double friction_coefficient;     // Friction coefficient
    double restitution_coefficient;  // Restitution coefficient
    // Optional: compute linear damping from restitution
    bool linear_use_restitution;     // If true, compute c from e
    double linear_en;                // Normal restitution for linear model
    double linear_et;                // Tangential restitution for linear model (-1 => use linear_en)
    
    // Rolling friction parameters
    double rolling_friction_coefficient;  // Rolling friction coefficient
    double rolling_stiffness;            // Rolling spring constant
    double rolling_damping;              // Rolling damping constant
    double rolling_torque_threshold;     // Threshold for rolling torque activation
    
    // Hertz model parameters
    double young_modulus;            // Young's modulus
    double poisson_ratio;            // Poisson's ratio
    
    // Sub-stepping parameters
    bool use_substeps;
    int max_substeps;
    
    // HYBRID: Adaptive collision detection thresholds
    int collision_simple_threshold;      // Triangle count below which to use sphere-sphere
    int collision_moderate_threshold;    // Triangle count below which to use triangle-triangle
    bool use_adaptive_collision;         // Enable adaptive algorithm selection
    
    // Contact history for tangential forces
    struct ContactHistory {
        Eigen::Vector3d tangential_overlap;
        bool in_contact;
        double last_update_time;
        // Enhanced contact history
        Eigen::Vector3d previous_normal;      // Previous contact normal
        Eigen::Vector3d tangential_spring;    // Tangential spring displacement
        double contact_duration;              // Duration of contact
    };
    map<pair<int, int>, ContactHistory> contact_history;
    
    // Clear contact history for pairs no longer in contact
    void update_contact_history(lexer *p);
    
    // Grid-based collision detection system
    sixdof_collision_grid* collision_grid;
    
    // NEW: Storage for collision forces and torques (for MPI communication)
    std::vector<Eigen::Vector3d> collision_forces;      // Collision forces for each object
    std::vector<Eigen::Vector3d> collision_torques;    // Collision torques for each object
    int max_objects;                                    // Maximum number of objects supported
    
    // NEW: AABB storage for fast pre-filtering
    std::vector<AABB> object_aabbs;                     // AABB for each object
    void update_aabbs(std::vector<sixdof_obj*> &fb_obj); // Update AABBs before collision detection
    
    // NEW: MPI communication functions
    void broadcast_collision_forces(lexer *p, ghostcell *pgc);
    void clear_collision_forces();
    void verify_collision_forces_synchronization(lexer *p, ghostcell *pgc); // Debug function
    
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
};

#endif 