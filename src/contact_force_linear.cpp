/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"contact_force_linear.h"
#include"contact_history.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include<cmath>

contact_force_linear::contact_force_linear(lexer *p)
{
    kn = p->R30;
    kt = p->R31;
    cn = 1.0e4;
    ct = 0.5e4;
    mu = p->R33;
    en = p->R32;
    et = -1.0;
    use_cor = true;
}

contact_force_linear::~contact_force_linear()
{
}

void contact_force_linear::compute(lexer *p, ghostcell*, sixdof_obj *obj1, sixdof_obj *obj2,
                                   const Eigen::Vector3d &contact_point,
                                   const Eigen::Vector3d &normal,
                                   double overlap,
                                   ContactHistory &history,
                                   Eigen::Vector3d &force, Eigen::Vector3d &torque)
{
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
    double v_rel_t_mag = v_rel_t.norm();

    double dt = p->simtime - history.t_last;
    if(dt <= 0.0) dt = p->dt;

    history.in_contact = true;
    history.t_last = p->simtime;

    double c_n = cn;
    if(use_cor)
    {
        const double en_lin = std::min(std::max(en, 1.0e-6), 0.999999);
        const double zeta = -std::log(en_lin) / std::sqrt(M_PI*M_PI + std::log(en_lin)*std::log(en_lin));
        const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
        c_n = 2.0 * zeta * std::sqrt(kn * meff);
    }

    double fn = kn * overlap - c_n * v_rel_n;
    fn = std::max(fn, 0.0);

    if(v_rel_t_mag > 1.0e-10)
    {
        history.s_t += v_rel_t * dt;
        history.s_t -= normal * normal.dot(history.s_t);

        double c_t = ct;
        if(use_cor)
        {
            double et_lin;
            if(et == -1.0)
            {
                const double en_lin = std::min(std::max(en, 1.0e-6), 0.999999);
                et_lin = en_lin;
            }
            else
            {
                et_lin = std::min(std::max(et, 1.0e-6), 0.999999);
            }
            const double zeta_t = -std::log(et_lin) / std::sqrt(M_PI*M_PI + std::log(et_lin)*std::log(et_lin));
            const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
            c_t = 2.0 * zeta_t * std::sqrt(kt * meff);
        }

        Eigen::Vector3d ft_vector = kt * history.s_t + c_t * v_rel_t;

        double ft_mag = ft_vector.norm();
        double ft_max = mu * fn;

        if(ft_mag > ft_max)
        {
            ft_vector *= (ft_max / ft_mag);
            history.s_t *= (ft_max / ft_mag);
        }

        force = fn * normal - ft_vector;
    }
    else
    {
        force = fn * normal;
    }

    torque = r1.cross(force);
}
