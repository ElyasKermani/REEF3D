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

#ifndef DEM_COLLISION_H_
#define DEM_COLLISION_H_

#include"contact_history.h"
#include"ground_contact.h"
#include"wall_contact.h"
#include<Eigen/Dense>
#include<vector>
#include<map>

class lexer;
class ghostcell;
class sixdof_obj;
class dem_collision_grid;
class contact_force;

using namespace std;

// Axis-aligned bounding box used for fast broad-phase rejection
struct AABB {
    Eigen::Vector3d min;
    Eigen::Vector3d max;

    AABB() : min(0,0,0), max(0,0,0) {}

    inline void update_from_sphere(const Eigen::Vector3d& center, double radius) {
        min = center - Eigen::Vector3d(radius, radius, radius);
        max = center + Eigen::Vector3d(radius, radius, radius);
    }

    inline bool overlaps(const AABB& other) const {
        return (max.x() >= other.min.x() && min.x() <= other.max.x()) &&
               (max.y() >= other.min.y() && min.y() <= other.max.y()) &&
               (max.z() >= other.min.z() && min.z() <= other.max.z());
    }

    // Inflate the box by a uniform margin
    inline void expand(double margin) {
        min.x() -= margin;
        min.y() -= margin;
        min.z() -= margin;
        max.x() += margin;
        max.y() += margin;
        max.z() += margin;
    }
};

// Available object-object contact-force models
enum class ContactForceModel {
    Linear,                       // Linear spring-dashpot
    Hertz,                        // Non-linear Hertzian elastic contact
    HertzMindlin,                 // Hertz with tangential history
    DemLinearLimited,             // linear spring-dashpot with tangential history rescaling at slip
    DemLinearNonLimited,          // linear spring-dashpot without tangential history rescaling at slip
    DemHertzianLimited,           // Hertz-Mindlin-style normal law; tangential history rescaling at slip
    DemHertzianNonLimited         // Hertz-Mindlin-style normal law; no tangential history rescaling at slip
};

enum class CollisionDetectionMode {
    SphereOnly,
    TriangleSATOnly,
    Adaptive
};

class dem_collision
{
public:

    dem_collision(lexer *p, ghostcell *pgc);
    virtual ~dem_collision();
    
    // Compute object-object collision forces and apply them as external loads
    void calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj);

    // Select the contact-force model used for object-object contacts
    void set_contact_force_model(ContactForceModel model);
    void set_adaptive_detection(bool enabled, int simple_tri, int moderate_tri);
    void set_detection_mode(CollisionDetectionMode mode) { detection_mode = mode; }
    
private:

    // Sphere-sphere narrow-phase test
    bool detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                         Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);

    // Triangle-mesh narrow-phase test using SAT (with optional BVH acceleration)
    bool detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                 Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);

    // Pick sphere-sphere or triangle-mesh based on object triangle counts
    bool detect_collision_adaptive(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                  Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap);

    // Triangle-triangle intersection test (separating-axis theorem)
    bool triangle_triangle_intersection(const Eigen::Vector3d v1[3], const Eigen::Vector3d v2[3],
                                      Eigen::Vector3d &contact, Eigen::Vector3d &normal, double &overlap);

    // Active contact-force kernel
    contact_force *p_force;
    lexer *p_lexer;

    ground_contact ground_;
    wall_contact   wall_;

    // Currently selected contact-force model
    ContactForceModel contact_model;

    // Rolling friction parameters
    double mu_r;        // rolling friction coefficient
    double kr;          // rolling spring stiffness
    double cr;          // rolling viscous damping
    double tau_r_max;   // rolling torque saturation threshold

    // Sub-stepping parameters
    bool use_substeps;
    int max_substeps;

    // Adaptive narrow-phase thresholds (triangle counts on each body)
    int simple_tri;       // use sphere-sphere if max(tri_i, tri_j) < this
    int moderate_tri;     // use triangle/SAT without BVH if below this; else BVH
    bool adaptive;       // if false, detection_mode alone selects narrow phase
    CollisionDetectionMode detection_mode;

    // Partner bounding-sphere radius multiplier for mesh-BVH pruning (from lexer R22)
    double bvh_prune_radius_scale;

    // Per-pair tangential history, keyed by (min(id1,id2), max(id1,id2))
    map<pair<int, int>, ContactHistory> contact_history;

    // Drop history entries for pairs that have been separated long enough
    void update_contact_history(lexer *p);

    // Spatial-hash grid used for broad-phase pair selection
    dem_collision_grid* collision_grid;

    // Per-object force/torque accumulators (computed on rank 0, broadcast to all)
    std::vector<Eigen::Vector3d> f_col;  // object–object contact force (world frame)
    std::vector<Eigen::Vector3d> t_col;  // object–object contact torque (world frame)
    int nobj;                            // number of floating bodies (from lexer X20 after inference)

    // AABBs for the currently active objects
    std::vector<AABB> aabbs;             // broad-phase AABB per object index
    void update_aabbs(std::vector<sixdof_obj*> &fb_obj);

    void broadcast_forces(lexer *p, ghostcell *pgc);
    void clear_forces();
    // Debug-only check that all ranks hold identical force/torque values
    void verify_sync(lexer *p, ghostcell *pgc);

    double calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2);

    // Resolve large overlaps by integrating the contact load over sub-time-steps
    void resolve_collision_with_substeps(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                       const Eigen::Vector3d &contact_point,
                                       const Eigen::Vector3d &normal,
                                       const double overlap,
                                       Eigen::Vector3d &force,
                                       Eigen::Vector3d &torque);

    // Symplectic velocity-Verlet update used inside the sub-step loop
    void velocity_verlet_step(lexer *p, ghostcell *pgc, sixdof_obj *obj,
                            const Eigen::Vector3d &force,
                            const Eigen::Vector3d &torque,
                            double dt);

    // Rolling friction torque opposing the relative angular velocity in the contact plane
    void calculate_rolling_friction_torque(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                         const Eigen::Vector3d &contact_point,
                                         const Eigen::Vector3d &normal,
                                         const double overlap,
                                         Eigen::Vector3d &rolling_torque);

    // Twisting resistance opposing relative angular velocity about the contact normal
    void calculate_twisting_resistance(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                     const Eigen::Vector3d &contact_point,
                                     const Eigen::Vector3d &normal,
                                     const double overlap,
                                     Eigen::Vector3d &twisting_torque);
};

#endif
