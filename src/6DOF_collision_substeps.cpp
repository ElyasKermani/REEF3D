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
#include<Eigen/Dense>
#include<iostream>

void sixdof_collision::resolve_collision_with_substeps(lexer *p, ghostcell *pgc, sixdof_obj *obj1, sixdof_obj *obj2,
                                                    const Eigen::Vector3d &contact_point, 
                                                    const Eigen::Vector3d &normal, 
                                                    const double overlap,
                                                    Eigen::Vector3d &force, 
                                                    Eigen::Vector3d &torque)
{
    // Number of sub-steps grows with the overlap severity, capped at max_substeps
    int num_substeps = std::min(static_cast<int>(ceil(overlap / (0.01 * std::min(obj1->radius, obj2->radius)))), max_substeps);

    if(p->mpirank==0 && p->count%p->P12==0)
    {
        cout<<"Using "<<num_substeps<<" sub-steps for collision resolution"<<endl;
    }

    double dt_sub = p->dt / static_cast<double>(num_substeps);

    // Snapshot the body states so we can roll them back at the end. The actual
    // force/torque is applied at the main time-step level by the caller; this
    // routine only previews the response over a finer sub-stepping schedule.
    Eigen::Vector3d orig_p1 = obj1->p_;
    Eigen::Vector3d orig_c1 = obj1->c_;
    Eigen::Vector3d orig_h1 = obj1->h_;
    Eigen::Vector4d orig_e1 = obj1->e_;

    Eigen::Vector3d orig_p2 = obj2->p_;
    Eigen::Vector3d orig_c2 = obj2->c_;
    Eigen::Vector3d orig_h2 = obj2->h_;
    Eigen::Vector4d orig_e2 = obj2->e_;

    for(int step = 0; step < num_substeps; ++step)
    {
        Eigen::Vector3d sub_force  = force  / static_cast<double>(num_substeps);
        Eigen::Vector3d sub_torque = torque / static_cast<double>(num_substeps);

        velocity_verlet_step(p, pgc, obj1, -sub_force, -sub_torque, dt_sub);
        velocity_verlet_step(p, pgc, obj2,  sub_force,  sub_torque, dt_sub);
    }

    // Restore the original states; the caller stores the load once via the
    // standard force accumulators to avoid double-counting.
    obj1->p_ = orig_p1;
    obj1->c_ = orig_c1;
    obj1->h_ = orig_h1;
    obj1->e_ = orig_e1;

    obj2->p_ = orig_p2;
    obj2->c_ = orig_c2;
    obj2->h_ = orig_h2;
    obj2->e_ = orig_e2;
}

void sixdof_collision::velocity_verlet_step(lexer *p, ghostcell *pgc, sixdof_obj *obj, 
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