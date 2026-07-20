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
#include<Eigen/Dense>
#include<iostream>

void dem_collision::snapshot_body_states(const vector<sixdof_obj*> &fb_obj,
                                         vector<BodySnapshotState> &snaps) const
{
    snaps.resize(fb_obj.size());

    for(size_t i = 0; i < fb_obj.size(); ++i)
    {
        snaps[i].p = fb_obj[i]->p_;
        snaps[i].c = fb_obj[i]->c_;
        snaps[i].h = fb_obj[i]->h_;
        snaps[i].e = fb_obj[i]->e_;
        snaps[i].omega_B = fb_obj[i]->omega_B;
        snaps[i].omega_I = fb_obj[i]->omega_I;
    }
}

void dem_collision::restore_body_states(lexer *p, vector<sixdof_obj*> &fb_obj,
                                        const vector<BodySnapshotState> &snaps)
{
    for(size_t i = 0; i < fb_obj.size(); ++i)
    {
        fb_obj[i]->p_ = snaps[i].p;
        fb_obj[i]->c_ = snaps[i].c;
        fb_obj[i]->h_ = snaps[i].h;
        fb_obj[i]->e_ = snaps[i].e;
        fb_obj[i]->omega_B = snaps[i].omega_B;
        fb_obj[i]->omega_I = snaps[i].omega_I;

        fb_obj[i]->quat_matrices(p);
    }
}

void dem_collision::integrate_contact_net_forces(lexer *p, ghostcell *pgc,
                                                 vector<sixdof_obj*> &fb_obj,
                                                 double dt_sub)
{
    for(size_t i = 0; i < fb_obj.size(); ++i)
    {
        velocity_verlet_step(p, pgc, fb_obj[i], f_col[i], t_col[i], dt_sub);
    }
}

void dem_collision::velocity_verlet_step(lexer *p, ghostcell *pgc, sixdof_obj *obj,
                                          const Eigen::Vector3d &force,
                                          const Eigen::Vector3d &torque,
                                          double dt)
{
    // Velocity-Verlet symplectic integrator: half kick, drift, half kick.
    // Linear momentum p_, angular momentum h_, position c_ and quaternion e_
    // are updated in-place.
    Eigen::Vector3d p_init = obj->p_;
    Eigen::Vector3d c_init = obj->c_;
    Eigen::Vector3d h_init = obj->h_;
    Eigen::Vector4d e_init = obj->e_;

    // First half-kick (linear and angular momentum)
    Eigen::Vector3d p_half = p_init + 0.5 * dt * force;
    Eigen::Vector3d torque_body = obj->Rinv_ * torque;
    Eigen::Vector3d h_half = h_init + 0.5 * dt * torque_body;

    // Drift: position from half-step linear momentum
    Eigen::Vector3d c_new = c_init + dt * p_half / obj->Mass_fb;

    // Quaternion update from body-frame angular velocity, dq/dt = 0.5 * G^T * omega_B
    Eigen::Vector3d omega_B = obj->I_.inverse() * h_half;
    obj->quat_matrices(p);
    Eigen::Vector4d de_dt = 0.5 * obj->G_.transpose() * omega_B;
    Eigen::Vector4d e_new = e_init + dt * de_dt;

    double e_norm = e_new.norm();
    if(e_norm > 1.0e-10)
        e_new /= e_norm;

    obj->e_ = e_new;
    obj->quat_matrices(p);

    // Second half-kick using the rotated body frame
    Eigen::Vector3d p_new = p_half + 0.5 * dt * force;
    torque_body = obj->Rinv_ * torque;
    Eigen::Vector3d h_new = h_half + 0.5 * dt * torque_body;

    obj->p_ = p_new;
    obj->c_ = c_new;
    obj->h_ = h_new;
    obj->e_ = e_new;

    obj->omega_B = obj->I_.inverse() * h_new;
    obj->omega_I = obj->R_ * obj->omega_B;
}
