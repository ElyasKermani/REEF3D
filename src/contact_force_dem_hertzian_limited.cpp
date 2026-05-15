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

#include"contact_force_dem_hertzian_limited.h"
#include"contact_history.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include<cmath>

contact_force_dem_hertzian_limited::contact_force_dem_hertzian_limited(lexer *p)
{
    E = p->R34;
    nu = p->R35;
    mu = p->R33;
    cor = p->R32;
}

contact_force_dem_hertzian_limited::~contact_force_dem_hertzian_limited()
{
}

void contact_force_dem_hertzian_limited::compute(lexer *p, ghostcell*, sixdof_obj *obj1, sixdof_obj *obj2,
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

    double dt = p->simtime - history.t_last;
    if(dt <= 0.0) dt = p->dt;

    history.in_contact = true;
    history.t_last = p->simtime;

    double R_eff = (obj1->radius * obj2->radius) / (obj1->radius + obj2->radius);
    double E_eff = effective_young_modulus(E, E, nu, nu);
    double G = E / (1.0 + nu);

    double mi = obj1->Mass_fb;
    double mj = obj2->Mass_fb;
    double meff = (mi*mj)/(mi+mj);

    double K_hertz = (4.0/3.0) * E_eff * std::sqrt(R_eff);
    const double en = std::min(std::max(cor, 1.0e-6), 0.999999);
    double ethan = -2.2664 * std::log(en) / std::sqrt(M_PI*M_PI + std::log(en)*std::log(en));

    Eigen::Vector3d Fn = (-(4.0/3.0) * E_eff * std::sqrt(R_eff) * std::pow(overlap, 1.5)
                          - std::sqrt(meff * K_hertz) * ethan * std::pow(overlap, 0.25) * v_rel_n) * normal;

    Eigen::Vector3d &overlap_t = history.s_t;

    if(v_rel_t.norm() > 1.0e-10)
    {
        overlap_t += v_rel_t * dt;
        overlap_t -= normal * normal.dot(overlap_t);

        double kt = 8.0 * G * std::sqrt(R_eff * overlap);
        Eigen::Vector3d Ft = -kt * overlap_t;
        double ft = Ft.norm();
        double ft_max = mu * Fn.norm();

        if(ft > ft_max)
        {
            if(overlap_t.norm() > 0.0)
            {
                Ft *= (ft_max/ft);
                overlap_t = -Ft / kt;
            }
            else
            {
                Ft.setZero();
            }
        }

        force = Fn + Ft;
    }
    else
    {
        force = Fn;
    }

    torque = r1.cross(force);
}
