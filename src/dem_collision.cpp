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

#include"dem_collision.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"dem_collision_grid.h"
#include"contact_force.h"
#include"contact_force_linear.h"
#include"contact_force_hertz.h"
#include"contact_force_hertz_mindlin.h"
#include"contact_force_dem_linear_limited.h"
#include"contact_force_dem_linear_nonlimited.h"
#include"contact_force_dem_hertzian_limited.h"
#include"contact_force_dem_hertzian_nonlimited.h"
#include<math.h>
#include<algorithm>
#include<iostream>

dem_collision::dem_collision(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"DEM collision model startup..."<<endl;
    
    p_lexer = p;
    p_force = 0;

    nobj = p->X20;
    f_col.resize(nobj);
    t_col.resize(nobj);
    aabbs.resize(nobj);
    clear_forces();

    // Default contact-force model
    set_contact_force_model(ContactForceModel::DemHertzianNonLimited);

    // Rolling/twisting friction parameters (lexer R40–R43)
    mu_r      = p->R40;
    kr        = p->R41;
    cr        = p->R42;
    tau_r_max = p->R43;

    use_substeps = (p->R51 != 0);
    max_substeps = std::max(1, p->R50);

    bvh_prune_radius_scale = std::max(p->R22, 1.0e-12);

    // Triangle-count thresholds for the adaptive narrow-phase test
    adaptive       = true;
    detection_mode = CollisionDetectionMode::Adaptive;
    simple_tri     = 50;    // below: sphere-sphere
    moderate_tri   = 500;   // below: triangle-mesh; above: triangle-mesh with BVH

    if(p->mpirank==0)
    {
        cout<<"DEM collision: Adaptive algorithm selection enabled"<<endl;
        cout<<"  Simple threshold (sphere-sphere): < "<<simple_tri<<" triangles"<<endl;
        cout<<"  Moderate threshold (triangle-mesh): < "<<moderate_tri<<" triangles"<<endl;
    }

    collision_grid = new dem_collision_grid(p, pgc);
}

dem_collision::~dem_collision()
{
    if(p_force)
        delete p_force;
    
    if(collision_grid)
        delete collision_grid;
}

void dem_collision::set_contact_force_model(ContactForceModel model)
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

    if(model==ContactForceModel::DemLinearLimited)
    p_force = new contact_force_dem_linear_limited(p_lexer);

    if(model==ContactForceModel::DemLinearNonLimited)
    p_force = new contact_force_dem_linear_nonlimited(p_lexer);

    if(model==ContactForceModel::DemHertzianLimited)
    p_force = new contact_force_dem_hertzian_limited(p_lexer);

    if(model==ContactForceModel::DemHertzianNonLimited)
    p_force = new contact_force_dem_hertzian_nonlimited(p_lexer);

    // Fall back to the linear model for any unrecognized enum value
    if(p_force==0)
    p_force = new contact_force_linear(p_lexer);

    if(p_lexer && p_lexer->mpirank==0)
    {
        cout<<"DEM collision: contact-force model = ";
        if(model==ContactForceModel::Linear)                            cout<<"Linear";
        if(model==ContactForceModel::Hertz)                             cout<<"Hertz";
        if(model==ContactForceModel::HertzMindlin)                      cout<<"Hertz-Mindlin";
        if(model==ContactForceModel::DemLinearLimited)           cout<<"DEM-Linear-Limited";
        if(model==ContactForceModel::DemLinearNonLimited)        cout<<"DEM-Linear-NonLimited";
        if(model==ContactForceModel::DemHertzianLimited)        cout<<"DEM-Hertzian-Limited";
        if(model==ContactForceModel::DemHertzianNonLimited)     cout<<"DEM-Hertzian-NonLimited";
        cout<<endl;
    }
}

void dem_collision::set_adaptive_detection(bool enabled, int simple_tri_in, int moderate_tri_in)
{
    adaptive = enabled;
    simple_tri = std::max(1, simple_tri_in);
    moderate_tri = std::max(simple_tri, moderate_tri_in);
}

void dem_collision::update_aabbs(vector<sixdof_obj*> &fb_obj)
{
    for(int i = 0; i < nobj; ++i)
    {
        aabbs[i].update_from_sphere(fb_obj[i]->c_, fb_obj[i]->radius);
    }
}

void dem_collision::calculate_collision_forces(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj,
                                               double dt_contact, bool finalize)
{
    nobj = static_cast<int>(fb_obj.size());
    if(nobj < 2)
    return;

    clear_forces();

    const int n_sub = estimate_substep_count(p, fb_obj, dt_contact);

    if(p->mpirank == 0)
    {
        update_contact_history(p);

        if(n_sub <= 1)
        {
            accumulate_object_contacts(p, pgc, fb_obj, dt_contact, finalize);
        }
        else
        {
            vector<BodySnapshotState> snaps;
            snapshot_body_states(fb_obj, snaps);

            vector<Eigen::Vector3d> f_sum(nobj, Eigen::Vector3d::Zero());
            vector<Eigen::Vector3d> t_sum(nobj, Eigen::Vector3d::Zero());

            const double dt_sub = dt_contact / static_cast<double>(n_sub);

            if(p->count%p->P12==0)
            {
                cout<<"DEM contact subcycling: "<<n_sub<<" steps, dt_sub="<<dt_sub
                    <<" s (dt_contact="<<dt_contact<<" s)"<<endl;
            }

            for(int sub = 0; sub < n_sub; ++sub)
            {
                const bool step_finalize = finalize && (sub == n_sub - 1);

                accumulate_object_contacts(p, pgc, fb_obj, dt_sub, step_finalize);

                for(int i = 0; i < nobj; ++i)
                {
                    f_sum[i] += f_col[i];
                    t_sum[i] += t_col[i];
                }

                integrate_contact_net_forces(p, pgc, fb_obj, dt_sub);
            }

            restore_body_states(p, fb_obj, snaps);

            const double inv_n = 1.0 / static_cast<double>(n_sub);
            for(int i = 0; i < nobj; ++i)
            {
                f_col[i] = f_sum[i] * inv_n;
                t_col[i] = t_sum[i] * inv_n;
            }
        }
    }

    broadcast_forces(p, pgc);

    // Apply the synchronized force/torque to each body on every rank
    for(int i = 0; i < nobj; ++i)
    {
        fb_obj[i]->Xext_dem = f_col[i](0);
        fb_obj[i]->Yext_dem = f_col[i](1);
        fb_obj[i]->Zext_dem = f_col[i](2);

        fb_obj[i]->Kext_dem = t_col[i](0);
        fb_obj[i]->Mext_dem = t_col[i](1);
        fb_obj[i]->Next_dem = t_col[i](2);
    }

    // Domain boundary contact (ground, walls) folds into Xext_dem/Yext_dem/Zext_dem so
    // hydrodynamic force routines that reset Xe/Ye/Ze do not erase these loads before solve_eqmotion.
    ground_.apply(p, pgc, fb_obj);
    wall_.apply(p, pgc, fb_obj);

    if(p->count%p->P38==0)
    {
        if(p->mpirank == 0)
        {
            cout<<"DEM collision: Applied forces to all objects on rank 0:"<<endl;
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
            cout<<"DEM collision: Rank "<<p->mpirank<<" received forces for "<<nobj<<" objects"<<endl;
        }

        // verify_sync(p, pgc); // Enable to check MPI sync
    }
}

int dem_collision::estimate_substep_count(lexer *p, const vector<sixdof_obj*> &fb_obj, double dt_contact) const
{
    if(!use_substeps || dt_contact <= 0.0 || fb_obj.size() < 2)
    return 1;

    const double kn = std::max(p->R30, 1.0);

    double m_min = fb_obj[0]->Mass_fb;
    for(size_t i = 0; i < fb_obj.size(); ++i)
    m_min = std::min(m_min, fb_obj[i]->Mass_fb);

    if(m_min <= 0.0)
    return 1;

    // Critical contact step ~ 0.2*sqrt(m/kn) (Cundall & Strack stability criterion)
    constexpr double beta = 0.2;
    const double dt_crit = beta * std::sqrt(m_min / kn);

    int n_stiff = 1;
    if(dt_crit > 1.0e-18)
    n_stiff = static_cast<int>(std::ceil(dt_contact / dt_crit));

    return std::max(1, std::min(n_stiff, max_substeps));
}

void dem_collision::accumulate_object_contacts(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj,
                                               double dt_step, bool finalize)
{
    clear_forces();

    update_aabbs(fb_obj);
    collision_grid->update_grid(p, pgc, fb_obj);

    const vector<pair<int, int>> potential_collisions =
        collision_grid->find_potential_collisions(p, pgc, fb_obj);

    if(p->count%p->P12==0 && potential_collisions.size() > 0)
    {
        cout<<"DEM collision: Found "<<potential_collisions.size()<<" potential collision pairs"<<endl;
    }

    for(const auto& pair : potential_collisions)
    {
        const int i = pair.first;
        const int j = pair.second;

        Eigen::Vector3d contact_point, normal;
        double overlap = 0.0;

        if(!aabbs[i].overlaps(aabbs[j]))
        continue;

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
        else if(adaptive)
        {
            collision_detected = detect_collision_adaptive(p, pgc, fb_obj[i], fb_obj[j],
                                                          contact_point, normal, overlap);
        }
        else
        {
            collision_detected = detect_collision(p, pgc, fb_obj[i], fb_obj[j],
                                                 contact_point, normal, overlap);
        }

        if(!collision_detected)
        continue;

        Eigen::Vector3d force, torque, rolling_torque, twisting_torque;

        const int id1 = fb_obj[i]->n6DOF;
        const int id2 = fb_obj[j]->n6DOF;
        const std::pair<int,int> key = std::make_pair(std::min(id1,id2), std::max(id1,id2));

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
                         dt_step, finalize,
                         force, torque);

        calculate_rolling_friction_torque(p, pgc, fb_obj[i], fb_obj[j],
                                        contact_point, normal, overlap,
                                        rolling_torque);

        calculate_twisting_resistance(p, pgc, fb_obj[i], fb_obj[j],
                                    contact_point, normal, overlap,
                                    twisting_torque);

        const Eigen::Vector3d r1 = contact_point - fb_obj[i]->c_;
        const Eigen::Vector3d r2 = contact_point - fb_obj[j]->c_;

        f_col[i]  -= force;
        f_col[j]  += force;

        t_col[i] += r1.cross(-force) - rolling_torque - twisting_torque;
        t_col[j] += r2.cross( force) + rolling_torque + twisting_torque;

        if(p->count%p->P12==0)
        {
            cout<<"DEM collision detected between objects "<<i<<" and "<<j<<endl;
            cout<<"  Model: ";
            switch(contact_model) {
                case ContactForceModel::Linear: cout<<"Linear"; break;
                case ContactForceModel::DemLinearLimited: cout<<"DEM-Linear-Limited"; break;
                case ContactForceModel::DemLinearNonLimited: cout<<"DEM-Linear-NonLimited"; break;
                case ContactForceModel::Hertz: cout<<"Hertz"; break;
                case ContactForceModel::HertzMindlin: cout<<"Hertz-Mindlin"; break;
                case ContactForceModel::DemHertzianLimited: cout<<"DEM-Hertzian-Limited"; break;
                case ContactForceModel::DemHertzianNonLimited: cout<<"DEM-Hertzian-NonLimited"; break;
                default: cout<<"Unknown"; break;
            }
            cout<<endl;
            cout<<"  Overlap: "<<overlap<<endl;
            cout<<"  Force: ["<<force(0)<<", "<<force(1)<<", "<<force(2)<<"]"<<endl;
        }
    }
}

bool dem_collision::detect_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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


void dem_collision::update_contact_history(lexer *p)
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

double dem_collision::calculate_distance_between_objects(sixdof_obj *obj1, sixdof_obj *obj2)
{
    return (obj2->c_ - obj1->c_).norm();
}

bool dem_collision::detect_triangle_collision(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                               Eigen::Vector3d &contact_point, Eigen::Vector3d &normal, double &overlap)
{
    // Bounding-sphere reject before paying the per-triangle cost
    if(!detect_collision(p, pgc, obj1, obj2, contact_point, normal, overlap))
        return false;

    // Optional BVH pruning: mesh BVH is in body frame (tri_x0). Transform the other
    // body's bounding sphere into local coordinates so the hierarchy stays valid as R, c change.
    if(obj1->use_bvh && obj1->mesh_bvh && obj1->mesh_bvh->is_built())
    {
        const Eigen::Vector3d c_loc = obj1->R_.transpose() * (obj2->c_ - obj1->c_);
        if(!obj1->mesh_bvh->intersects_sphere(c_loc, obj2->radius * bvh_prune_radius_scale))
        {
            return false;
        }
    }

    if(obj2->use_bvh && obj2->mesh_bvh && obj2->mesh_bvh->is_built())
    {
        const Eigen::Vector3d c_loc = obj2->R_.transpose() * (obj1->c_ - obj2->c_);
        if(!obj2->mesh_bvh->intersects_sphere(c_loc, obj1->radius * bvh_prune_radius_scale))
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

bool dem_collision::triangle_triangle_intersection(const Eigen::Vector3d v1[3], const Eigen::Vector3d v2[3],
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

bool dem_collision::detect_collision_adaptive(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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

void dem_collision::calculate_rolling_friction_torque(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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

void dem_collision::calculate_twisting_resistance(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
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

void dem_collision::clear_forces()
{
    for(int i = 0; i < nobj; ++i)
    {
        f_col[i].setZero();
        t_col[i].setZero();
    }
}

void dem_collision::broadcast_forces(lexer *p, ghostcell *pgc)
{
    for(int i = 0; i < nobj; ++i)
    {
        MPI_Bcast(&f_col[i](0),  3, MPI_DOUBLE, 0, pgc->mpi_comm);
        MPI_Bcast(&t_col[i](0), 3, MPI_DOUBLE, 0, pgc->mpi_comm);
    }

    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"DEM collision: Broadcasted forces and torques to all processors"<<endl;
    }
}

// Diagnostic helper: confirm that every rank holds the same force/torque values.
// Not part of the normal force application path; call from a debug session if
// MPI synchronization is suspected.
void dem_collision::verify_sync(lexer *p, ghostcell *pgc)
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
            cout<<"DEM collision: SUCCESS - All collision forces are synchronized across processors"<<endl;
        else
            cout<<"DEM collision: ERROR - Collision forces are NOT synchronized across processors!"<<endl;
    }
}