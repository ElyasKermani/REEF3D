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
#include"contact_force.h"
#include"contact_force_linear.h"
#include"contact_force_hertz.h"
#include"contact_force_hertz_mindlin.h"
#include"contact_force_phasicflow_linear_limited.h"
#include"contact_force_phasicflow_linear_nonlimited.h"
#include"contact_force_phasicflow_nonlinear_limited.h"
#include"contact_force_phasicflow_nonlinear_nonlimited.h"
#include<math.h>
#include<iostream>

sixdof_collision::sixdof_collision(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF Collision Model startup..."<<endl;
    
    p_lexer = p;
    p_force = 0;

    nobj = p->X20;
    f_col.resize(nobj);
    t_col.resize(nobj);
    aabbs.resize(nobj);
    clear_forces();

    // Default contact-force model
    set_contact_force_model(ContactForceModel::PhasicFlowNonLinearNonLimited);

    // Rolling/twisting friction parameters
    mu_r      = 0.2;
    kr        = 0.5e6;   // [N m / rad]
    cr        = 0.5e4;   // [N m s / rad]
    tau_r_max = 1.0e-3;  // [N m]

    use_substeps = true;
    max_substeps = 10;

    // Triangle-count thresholds for the adaptive narrow-phase test
    adaptive       = true;
    detection_mode = CollisionDetectionMode::Adaptive;
    simple_tri     = 50;    // below: sphere-sphere
    moderate_tri   = 500;   // below: triangle-mesh; above: triangle-mesh with BVH

    if(p->mpirank==0)
    {
        cout<<"6DOF Collision: Adaptive algorithm selection enabled"<<endl;
        cout<<"  Simple threshold (sphere-sphere): < "<<simple_tri<<" triangles"<<endl;
        cout<<"  Moderate threshold (triangle-mesh): < "<<moderate_tri<<" triangles"<<endl;
    }

    collision_grid = new sixdof_collision_grid(p, pgc);
}

sixdof_collision::~sixdof_collision()
{
    if(p_force)
        delete p_force;
    
    if(collision_grid)
        delete collision_grid;
}

void sixdof_collision::set_contact_force_model(ContactForceModel model)
{
    // Skip re-allocation when the requested model is already active
    if(p_force!=0 && contact_model==model)
    return;

    contact_model = model;

    if(p_force)
    {
        delete p_force;
        p_force = 0;
    }

    if(model==ContactForceModel::Linear)
    p_force = new contact_force_linear(p_lexer);

    if(model==ContactForceModel::Hertz)
    p_force = new contact_force_hertz(p_lexer);

    if(model==ContactForceModel::HertzMindlin)
    p_force = new contact_force_hertz_mindlin(p_lexer);

    if(model==ContactForceModel::PhasicFlowLinearLimited)
    p_force = new contact_force_phasicflow_linear_limited(p_lexer);

    if(model==ContactForceModel::PhasicFlowLinearNonLimited)
    p_force = new contact_force_phasicflow_linear_nonlimited(p_lexer);

    if(model==ContactForceModel::PhasicFlowNonLinearLimited)
    p_force = new contact_force_phasicflow_nonlinear_limited(p_lexer);

    if(model==ContactForceModel::PhasicFlowNonLinearNonLimited)
    p_force = new contact_force_phasicflow_nonlinear_nonlimited(p_lexer);

    // Fall back to the linear model for any unrecognized enum value
    if(p_force==0)
    p_force = new contact_force_linear(p_lexer);

    if(p_lexer && p_lexer->mpirank==0)
    {
        cout<<"6DOF Collision: contact-force model = ";
        if(model==ContactForceModel::Linear)                            cout<<"Linear";
        if(model==ContactForceModel::Hertz)                             cout<<"Hertz";
        if(model==ContactForceModel::HertzMindlin)                      cout<<"Hertz-Mindlin";
        if(model==ContactForceModel::PhasicFlowLinearLimited)           cout<<"PhasicFlow-Linear-Limited";
        if(model==ContactForceModel::PhasicFlowLinearNonLimited)        cout<<"PhasicFlow-Linear-NonLimited";
        if(model==ContactForceModel::PhasicFlowNonLinearLimited)        cout<<"PhasicFlow-NonLinear-Limited";
        if(model==ContactForceModel::PhasicFlowNonLinearNonLimited)     cout<<"PhasicFlow-NonLinear-NonLimited";
        cout<<endl;
    }
}

void sixdof_collision::set_adaptive_detection(bool enabled, int simple_tri_in, int moderate_tri_in)
{
    adaptive = enabled;
    simple_tri = std::max(1, simple_tri_in);
    moderate_tri = std::max(simple_tri, moderate_tri_in);
}

void sixdof_collision::update_aabbs(vector<sixdof_obj*> &fb_obj)
{
    for(int i = 0; i < nobj; ++i)
    {
        aabbs[i].update_from_sphere(fb_obj[i]->c_, fb_obj[i]->radius);
    }
}

void sixdof_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    clear_forces();

    // Rank 0 owns the contact-force computation; results are broadcast below
    if(p->mpirank == 0)
    {
        update_aabbs(fb_obj);
        update_contact_history(p);
        collision_grid->update_grid(p, pgc, fb_obj);

        std::vector<std::pair<int, int>> potential_collisions =
            collision_grid->find_potential_collisions(p, pgc, fb_obj);

        if(p->count%p->P12==0 && potential_collisions.size() > 0)
        {
            cout<<"6DOF Collision: Found "<<potential_collisions.size()<<" potential collision pairs"<<endl;
        }

        for(const auto& pair : potential_collisions)
        {
            int i = pair.first;
            int j = pair.second;

            Eigen::Vector3d contact_point, normal;
            double overlap = 0.0;

            // Cheap AABB rejection before the narrow-phase test
            if(!aabbs[i].overlaps(aabbs[j]))
            {
                continue;
            }

            bool collision_detected = false;
            
            if(detection_mode == CollisionDetectionMode::SphereOnly)
            {
                collision_detected = detect_collision(p, pgc, fb_obj[i], fb_obj[j],
                                                     contact_point, normal, overlap);
            }
            else if(detection_mode == CollisionDetectionMode::TriangleSATOnly)
            {
                collision_detected = detect_triangle_collision(p, pgc, fb_obj[i], fb_obj[j],
                                                              contact_point, normal, overlap);
            }
            else
            {
                if(adaptive)
                {
                    collision_detected = detect_collision_adaptive(p, pgc, fb_obj[i], fb_obj[j],
                                                                  contact_point, normal, overlap);
                }
                else
                {
                    collision_detected = detect_collision(p, pgc, fb_obj[i], fb_obj[j],
                                                         contact_point, normal, overlap);
                }
            }
            
            if(collision_detected)
            {
                Eigen::Vector3d force, torque, rolling_torque, twisting_torque;

                const int id1 = fb_obj[i]->n6DOF;
                const int id2 = fb_obj[j]->n6DOF;
                const std::pair<int,int> key = std::make_pair(std::min(id1,id2), std::max(id1,id2));

                // Initialize a fresh contact-history entry the first time this pair touches
                if(contact_history.find(key) == contact_history.end())
                {
                    ContactHistory h;
                    h.s_t.setZero();
                    h.in_contact = true;
                    h.t_last = p->simtime;
                    contact_history[key] = h;
                }

                p_force->compute(p, pgc, fb_obj[i], fb_obj[j],
                                 contact_point, normal, overlap,
                                 contact_history[key],
                                 force, torque);

                calculate_rolling_friction_torque(p, pgc, fb_obj[i], fb_obj[j],
                                                contact_point, normal, overlap,
                                                rolling_torque);

                calculate_twisting_resistance(p, pgc, fb_obj[i], fb_obj[j],
                                            contact_point, normal, overlap,
                                            twisting_torque);

                if(use_substeps && overlap > 0.01 * fb_obj[i]->radius)
                {
                    resolve_collision_with_substeps(p, pgc, fb_obj[i], fb_obj[j],
                                                 contact_point, normal, overlap,
                                                 force, torque);
                }

                // Apply force/torque pair: action +force on obj j, reaction -force on obj i.
                // Each contact torque uses its own lever arm; rolling/twisting torques are
                // returned in the obj j frame and flipped on obj i.
                const Eigen::Vector3d r1 = contact_point - fb_obj[i]->c_;
                const Eigen::Vector3d r2 = contact_point - fb_obj[j]->c_;

                f_col[i]  -= force;
                f_col[j]  += force;

                t_col[i] += r1.cross(-force) - rolling_torque - twisting_torque;
                t_col[j] += r2.cross( force) + rolling_torque + twisting_torque;

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

    broadcast_forces(p, pgc);

    // Apply the synchronized force/torque to each body on every rank
    for(int i = 0; i < nobj; ++i)
    {
        fb_obj[i]->Xext += f_col[i](0);
        fb_obj[i]->Yext += f_col[i](1);
        fb_obj[i]->Zext += f_col[i](2);

        fb_obj[i]->Kext += t_col[i](0);
        fb_obj[i]->Mext += t_col[i](1);
        fb_obj[i]->Next += t_col[i](2);
    }

    calculate_ground_contact_forces(p, pgc, fb_obj);
    calculate_boundary_wall_contact_forces(p, pgc, fb_obj);

    if(p->count%p->P38==0)
    {
        if(p->mpirank == 0)
        {
            cout<<"6DOF Collision: Applied forces to all objects on rank 0:"<<endl;
            for(int i = 0; i < nobj; ++i)
            {
                if(f_col[i].norm() > 1.0e-10 || t_col[i].norm() > 1.0e-10)
                {
                    cout<<"  Object "<<i<<": Force=["<<f_col[i](0)<<", "<<f_col[i](1)<<", "<<f_col[i](2)<<"]"<<endl;
                    cout<<"            Torque=["<<t_col[i](0)<<", "<<t_col[i](1)<<", "<<t_col[i](2)<<"]"<<endl;
                }
            }
        }
        else
        {
            cout<<"6DOF Collision: Rank "<<p->mpirank<<" received forces for "<<nobj<<" objects"<<endl;
        }

        // verify_sync(p, pgc); // Enable to check MPI sync
    }
}

bool sixdof_collision::detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                      Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Sphere-sphere intersection test using the bounding radii of each body
    Eigen::Vector3d center1 = obj1->c_;
    Eigen::Vector3d center2 = obj2->c_;

    Eigen::Vector3d center_diff = center2 - center1;
    double distance = center_diff.norm();
    double sum_radii = obj1->radius + obj2->radius;

    if(distance < sum_radii)
    {
        overlap = sum_radii - distance;

        if(distance > 1.0e-10)
        {
            normal = center_diff / distance;
        }
        else
        {
            // Fallback when the two centres coincide
            normal = Eigen::Vector3d(0.0, 0.0, 1.0);
        }

        // Contact point at the midpoint of the penetration along the normal
        contact_point = center1 + normal * (obj1->radius - 0.5 * overlap);

        return true;
    }

    return false;
}


void sixdof_collision::update_contact_history(lexer *p)
{
    // Drop pairs that have been separated for longer than the timeout, and reset
    // the in_contact flag on the remaining entries; it is set to true again by
    // the dispatch loop whenever a contact is detected this step.
    const double contact_timeout = 1.0;

    for(auto it = contact_history.begin(); it != contact_history.end();)
    {
        if(!it->second.in_contact && (p->simtime - it->second.t_last) > contact_timeout)
        {
            it = contact_history.erase(it);
        }
        else
        {
            it->second.in_contact = false;
            ++it;
        }
    }
}

double sixdof_collision::calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2)
{
    return (obj2->c_ - obj1->c_).norm();
}

bool sixdof_collision::detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                               Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Bounding-sphere reject before paying the per-triangle cost
    if(!detect_collision(p, pgc, obj1, obj2, contact_point, normal, overlap))
        return false;

    // Optional BVH reject: skip the full pairwise sweep if the other body's
    // bounding sphere does not intersect our triangle BVH
    if(obj1->use_bvh && obj1->mesh_bvh && obj1->mesh_bvh->is_built())
    {
        if(!obj1->mesh_bvh->intersects_sphere(obj2->c_, obj2->radius))
        {
            return false;
        }
    }

    if(obj2->use_bvh && obj2->mesh_bvh && obj2->mesh_bvh->is_built())
    {
        if(!obj2->mesh_bvh->intersects_sphere(obj1->c_, obj1->radius))
        {
            return false;
        }
    }

    double min_overlap = 1e10;
    bool collision_found = false;

    // Transform each triangle into world coordinates and test all pairs
    for(int i=0; i<obj1->tricount; ++i)
    {
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
            Eigen::Vector3d v2[3];
            for(int q=0; q<3; ++q)
            {
                Eigen::Vector3d local_point(obj2->tri_x0[j][q], obj2->tri_y0[j][q], obj2->tri_z0[j][q]);
                Eigen::Vector3d global_point;
                obj2->motionext_trans(p, pgc, local_point, global_point);
                v2[q] = global_point;
            }

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
    // Separating-axis theorem with 11 candidate axes:
    //   2 face normals + 9 edge-edge cross products.
    // Two convex shapes intersect iff none of these axes separate them.
    const double EPSILON = 1e-8;

    Eigen::Vector3d n1 = (v1[1] - v1[0]).cross(v1[2] - v1[0]);
    Eigen::Vector3d n2 = (v2[1] - v2[0]).cross(v2[2] - v2[0]);

    // Reject degenerate (zero-area) triangles
    double n1_len = n1.norm();
    double n2_len = n2.norm();
    if(n1_len < EPSILON || n2_len < EPSILON)
        return false;

    n1 /= n1_len;
    n2 /= n2_len;

    Eigen::Vector3d e1[3], e2[3];
    e1[0] = v1[1] - v1[0];
    e1[1] = v1[2] - v1[1];
    e1[2] = v1[0] - v1[2];

    e2[0] = v2[1] - v2[0];
    e2[1] = v2[2] - v2[1];
    e2[2] = v2[0] - v2[2];

    double min_penetration = 1e10;
    Eigen::Vector3d best_axis;
    bool found_separation = false;

    // Project both triangles onto a candidate axis; return false on separation.
    // Otherwise update the running minimum-penetration axis and continue.
    auto test_axis = [&](const Eigen::Vector3d& axis) -> bool
    {
        double axis_len = axis.norm();
        if(axis_len < EPSILON)
            return true;

        Eigen::Vector3d axis_norm = axis / axis_len;

        double proj1[3];
        for(int i = 0; i < 3; ++i)
            proj1[i] = v1[i].dot(axis_norm);

        double min1 = std::min({proj1[0], proj1[1], proj1[2]});
        double max1 = std::max({proj1[0], proj1[1], proj1[2]});

        double proj2[3];
        for(int i = 0; i < 3; ++i)
            proj2[i] = v2[i].dot(axis_norm);

        double min2 = std::min({proj2[0], proj2[1], proj2[2]});
        double max2 = std::max({proj2[0], proj2[1], proj2[2]});

        if(max1 < min2 - EPSILON || max2 < min1 - EPSILON)
        {
            found_separation = true;
            return false;
        }

        double penetration = std::min(max1 - min2, max2 - min1);

        if(penetration < min_penetration)
        {
            min_penetration = penetration;
            best_axis = axis_norm;

            // Orient the contact normal from triangle 2 toward triangle 1
            Eigen::Vector3d center1 = (v1[0] + v1[1] + v1[2]) / 3.0;
            Eigen::Vector3d center2 = (v2[0] + v2[1] + v2[2]) / 3.0;
            if((center1 - center2).dot(best_axis) < 0)
                best_axis = -best_axis;
        }

        return true;
    };

    // Face normals
    if(!test_axis(n1)) return false;
    if(!test_axis(n2)) return false;

    // Edge-edge cross products
    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            Eigen::Vector3d axis = e1[i].cross(e2[j]);
            if(!test_axis(axis)) return false;
        }
    }

    overlap = min_penetration;
    normal = best_axis;

    // Approximate the contact point as the centroid of vertices lying close to
    // the opposing triangle's plane; fall back to the centroid midpoint.
    std::vector<Eigen::Vector3d> contact_points;

    for(int i = 0; i < 3; ++i)
    {
        double dist = (v1[i] - v2[0]).dot(n2);
        if(fabs(dist) < min_penetration + EPSILON)
            contact_points.push_back(v1[i]);

        dist = (v2[i] - v1[0]).dot(n1);
        if(fabs(dist) < min_penetration + EPSILON)
            contact_points.push_back(v2[i]);
    }

    if(contact_points.empty())
    {
        Eigen::Vector3d center1 = (v1[0] + v1[1] + v1[2]) / 3.0;
        Eigen::Vector3d center2 = (v2[0] + v2[1] + v2[2]) / 3.0;
        contact = 0.5 * (center1 + center2);
    }
    else
    {
        contact = Eigen::Vector3d::Zero();
        for(const auto& pt : contact_points)
            contact += pt;
        contact /= contact_points.size();
    }

    return true;
}

bool sixdof_collision::detect_collision_adaptive(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                 Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Pick the narrow-phase algorithm based on the more complex of the two bodies.
    int max_tri_count = std::max(obj1->tricount, obj2->tricount);

    if(max_tri_count < simple_tri)
    {
        // Cheap bounding-sphere approximation for low-poly geometry
        return detect_collision(p, pgc, obj1, obj2, contact_point, normal, overlap);
    }

    // Otherwise fall back to the SAT triangle-triangle test (BVH-accelerated when available)
    return detect_triangle_collision(p, pgc, obj1, obj2, contact_point, normal, overlap);
}

void sixdof_collision::calculate_rolling_friction_torque(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                      const Eigen::Vector3d &contact_point,
                                                      const Eigen::Vector3d &normal,
                                                      const double overlap,
                                                      Eigen::Vector3d &rolling_torque)
{
    // Component of the relative angular velocity that lies in the contact plane
    Eigen::Vector3d omega_rel   = obj2->omega_I - obj1->omega_I;
    Eigen::Vector3d omega_rel_t = omega_rel - normal * normal.dot(omega_rel);

    double rolling_velocity = omega_rel_t.norm();

    if(rolling_velocity > 1.0e-10)
    {
        Eigen::Vector3d rolling_direction = omega_rel_t / rolling_velocity;

        double rolling_torque_mag = mu_r * overlap *
                                  (kr * rolling_velocity +
                                   cr   * rolling_velocity * rolling_velocity);

        if(rolling_torque_mag > tau_r_max)
            rolling_torque_mag = tau_r_max;

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
    // Component of the relative angular velocity along the contact normal
    Eigen::Vector3d omega_rel  = obj2->omega_I - obj1->omega_I;
    double twisting_velocity   = normal.dot(omega_rel);

    if(std::abs(twisting_velocity) > 1.0e-10)
    {
        double twisting_torque_mag = mu_r * overlap *
                                   (kr * std::abs(twisting_velocity) +
                                    cr   * twisting_velocity * twisting_velocity);

        if(twisting_torque_mag > tau_r_max)
            twisting_torque_mag = tau_r_max;

        twisting_torque = -twisting_torque_mag * normal * (twisting_velocity > 0 ? 1.0 : -1.0);
    }
    else
    {
        twisting_torque.setZero();
    }
}

void sixdof_collision::clear_forces()
{
    for(int i = 0; i < nobj; ++i)
    {
        f_col[i].setZero();
        t_col[i].setZero();
    }
}

void sixdof_collision::broadcast_forces(lexer *p, ghostcell *pgc)
{
    for(int i = 0; i < nobj; ++i)
    {
        MPI_Bcast(&f_col[i](0),  3, MPI_DOUBLE, 0, pgc->mpi_comm);
        MPI_Bcast(&t_col[i](0), 3, MPI_DOUBLE, 0, pgc->mpi_comm);
    }

    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"6DOF Collision: Broadcasted forces and torques to all processors"<<endl;
    }
}

// Diagnostic helper: confirm that every rank holds the same force/torque values.
// Not part of the normal force application path; call from a debug session if
// MPI synchronization is suspected.
void sixdof_collision::verify_sync(lexer *p, ghostcell *pgc)
{
    bool forces_synchronized = true;

    for(int i = 0; i < nobj; ++i)
    {
        double force_x_sum = 0.0, force_y_sum = 0.0, force_z_sum = 0.0;
        double torque_x_sum = 0.0, torque_y_sum = 0.0, torque_z_sum = 0.0;

        MPI_Allreduce(&f_col[i](0), &force_x_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&f_col[i](1), &force_y_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&f_col[i](2), &force_z_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);

        MPI_Allreduce(&t_col[i](0), &torque_x_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&t_col[i](1), &torque_y_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(&t_col[i](2), &torque_z_sum, 1, MPI_DOUBLE, MPI_SUM, pgc->mpi_comm);

        int num_procs;
        MPI_Comm_size(pgc->mpi_comm, &num_procs);

        double expected_force_x  = f_col[i](0)  * num_procs;
        double expected_force_y  = f_col[i](1)  * num_procs;
        double expected_force_z  = f_col[i](2)  * num_procs;

        double expected_torque_x = t_col[i](0) * num_procs;
        double expected_torque_y = t_col[i](1) * num_procs;
        double expected_torque_z = t_col[i](2) * num_procs;

        const double tolerance = 1.0e-10;
        if(fabs(force_x_sum  - expected_force_x ) > tolerance ||
           fabs(force_y_sum  - expected_force_y ) > tolerance ||
           fabs(force_z_sum  - expected_force_z ) > tolerance ||
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
            cout<<"6DOF Collision: SUCCESS - All collision forces are synchronized across processors"<<endl;
        else
            cout<<"6DOF Collision: ERROR - Collision forces are NOT synchronized across processors!"<<endl;
    }
}

// ============================================================================
// GROUND CONTACT FORCES
// ============================================================================

void sixdof_collision::calculate_ground_contact_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Spring-damper-friction model shared with the mooring system in
    // src/mooring_dynamic_rhs.cpp.
    const double zg = 0.0;        // Seabed elevation
    const double Kg = 3.0e6;      // Ground stiffness         [N / m^2]
    const double mu = 0.3;        // Coulomb friction coefficient
    const double xi = 1.0;        // Damping ratio (critical damping)
    const double nu = 0.01;       // Friction velocity scale  [m / s]

    int obj_count = fb_obj.size();

    for(int i = 0; i < obj_count; i++)
    {
        sixdof_obj* obj = fb_obj[i];

        // Conservative lowest-point estimate from the bounding radius
        double z_bottom = obj->c_(2) - obj->radius;

        if((zg - z_bottom) > 0.0)
        {
            double overlap = zg - z_bottom;

            // Normal force: spring proportional to overlap with one-way damping
            // (damping only acts while the body is moving into the ground).
            double k_eff   = Kg * obj->radius;
            double m_eff   = obj->Mass_fb;
            double c_eff   = 2.0 * xi * sqrt(k_eff * m_eff);
            double v_normal = obj->u_fb(2);
            double F_normal = k_eff * overlap - c_eff * max(v_normal, 0.0);

            obj->Ze += F_normal;

            // Tangential friction: smoothed Coulomb model that transitions
            // linearly through |v|<nu and saturates at mu * F_normal beyond.
            double v_x = obj->u_fb(0);
            double v_y = obj->u_fb(1);
            double v_horiz = sqrt(v_x*v_x + v_y*v_y);

            double vx_norm = v_x / max(nu, v_horiz);
            double vy_norm = v_y / max(nu, v_horiz);

            double F_friction_mag = mu * F_normal;
            double Fx_friction = -F_friction_mag * sin(PI/2.0 * vx_norm);
            double Fy_friction = -F_friction_mag * sin(PI/2.0 * vy_norm);

            obj->Xe += Fx_friction;
            obj->Ye += Fy_friction;

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Ground Contact - Object " << i
                     << ": z_bottom=" << z_bottom
                     << " overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << " F_friction=" << sqrt(Fx_friction*Fx_friction + Fy_friction*Fy_friction) << " N"
                     << " v_z=" << v_normal << " m/s"
                     << endl;
            }
        }
    }
}

// ============================================================================
// BOUNDARY WALL CONTACT FORCES (Side Walls)
// ============================================================================

void sixdof_collision::calculate_boundary_wall_contact_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj)
{
    // Same spring-damper-friction model as the seabed contact
    const double K_wall = 3.0e6;  // Wall stiffness        [N / m^2]
    const double mu     = 0.3;    // Coulomb friction coefficient
    const double xi     = 1.0;    // Damping ratio
    const double nu     = 0.01;   // Friction velocity scale [m / s]

    int obj_count = fb_obj.size();

    for(int i = 0; i < obj_count; i++)
    {
        sixdof_obj* obj = fb_obj[i];

        double x_center = obj->c_(0);
        double y_center = obj->c_(1);
        double radius   = obj->radius;

        // ----- X walls -----------------------------------------------------

        double x_left = x_center - radius;
        if(x_left < p->xcoormin)
        {
            double overlap = p->xcoormin - x_left;

            double k_eff   = K_wall * radius;
            double m_eff   = obj->Mass_fb;
            double c_eff   = 2.0 * xi * sqrt(k_eff * m_eff);
            double v_normal = obj->u_fb(0);
            double F_normal = k_eff * overlap - c_eff * min(v_normal, 0.0);

            obj->Xe += F_normal;

            double v_y = obj->u_fb(1);
            double v_z = obj->u_fb(2);
            double v_tangential = sqrt(v_y*v_y + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                double F_friction_mag = mu * F_normal;
                double vy_norm = v_y / max(nu, v_tangential);
                double vz_norm = v_z / max(nu, v_tangential);

                obj->Ye -= F_friction_mag * sin(PI/2.0 * vy_norm);
                obj->Ze -= F_friction_mag * sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (X-min) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        double x_right = x_center + radius;
        if(x_right > p->xcoormax)
        {
            double overlap = x_right - p->xcoormax;

            double k_eff   = K_wall * radius;
            double m_eff   = obj->Mass_fb;
            double c_eff   = 2.0 * xi * sqrt(k_eff * m_eff);
            double v_normal = obj->u_fb(0);
            double F_normal = k_eff * overlap - c_eff * max(v_normal, 0.0);

            obj->Xe -= F_normal;

            double v_y = obj->u_fb(1);
            double v_z = obj->u_fb(2);
            double v_tangential = sqrt(v_y*v_y + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                double F_friction_mag = mu * F_normal;
                double vy_norm = v_y / max(nu, v_tangential);
                double vz_norm = v_z / max(nu, v_tangential);

                obj->Ye -= F_friction_mag * sin(PI/2.0 * vy_norm);
                obj->Ze -= F_friction_mag * sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (X-max) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        // ----- Y walls -----------------------------------------------------

        double y_front = y_center - radius;
        if(y_front < p->ycoormin)
        {
            double overlap = p->ycoormin - y_front;

            double k_eff   = K_wall * radius;
            double m_eff   = obj->Mass_fb;
            double c_eff   = 2.0 * xi * sqrt(k_eff * m_eff);
            double v_normal = obj->u_fb(1);
            double F_normal = k_eff * overlap - c_eff * min(v_normal, 0.0);

            obj->Ye += F_normal;

            double v_x = obj->u_fb(0);
            double v_z = obj->u_fb(2);
            double v_tangential = sqrt(v_x*v_x + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                double F_friction_mag = mu * F_normal;
                double vx_norm = v_x / max(nu, v_tangential);
                double vz_norm = v_z / max(nu, v_tangential);

                obj->Xe -= F_friction_mag * sin(PI/2.0 * vx_norm);
                obj->Ze -= F_friction_mag * sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (Y-min) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        double y_back = y_center + radius;
        if(y_back > p->ycoormax)
        {
            double overlap = y_back - p->ycoormax;

            double k_eff   = K_wall * radius;
            double m_eff   = obj->Mass_fb;
            double c_eff   = 2.0 * xi * sqrt(k_eff * m_eff);
            double v_normal = obj->u_fb(1);
            double F_normal = k_eff * overlap - c_eff * max(v_normal, 0.0);

            obj->Ye -= F_normal;

            double v_x = obj->u_fb(0);
            double v_z = obj->u_fb(2);
            double v_tangential = sqrt(v_x*v_x + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                double F_friction_mag = mu * F_normal;
                double vx_norm = v_x / max(nu, v_tangential);
                double vz_norm = v_z / max(nu, v_tangential);

                obj->Xe -= F_friction_mag * sin(PI/2.0 * vx_norm);
                obj->Ze -= F_friction_mag * sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (Y-max) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }
    }
}