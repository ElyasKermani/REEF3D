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

#include"contact_force_hertz_mindlin.h"
#include"contact_history.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include<cmath>

contact_force_hertz_mindlin::contact_force_hertz_mindlin(lexer *p)
{
    E = p->R34;
    nu = p->R35;
    mu = p->R33;
    cor = p->R32;
}

contact_force_hertz_mindlin::~contact_force_hertz_mindlin()
{
}

void contact_force_hertz_mindlin::compute(lexer *p, ghostcell*, sixdof_obj *obj1, sixdof_obj *obj2,
                                          const Eigen::Vector3d &contact_point,
                                          const Eigen::Vector3d &normal,
                                          double overlap,
                                          ContactHistory &history,
                                          double dt_contact,
                                          bool finalize,
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

    history.in_contact = true;

    const double R_eff = effective_radius(obj1->radius, obj2->radius);
    const double E_eff = effective_young_modulus(E, E, nu, nu);

    const double k_hertz = hertz_stiffness(E_eff, R_eff);
    const double fn_elastic = (4.0/3.0) * k_hertz * std::pow(overlap, 1.5);
    const double Sn = 2.0 * E_eff * std::sqrt(R_eff * overlap);
    const double meff = (obj1->Mass_fb * obj2->Mass_fb) / (obj1->Mass_fb + obj2->Mass_fb);
    const double en_ = std::min(std::max(cor, 1.0e-6), 0.999999);
    const double beta = std::log(en_) / std::sqrt(M_PI*M_PI + std::log(en_)*std::log(en_));
    const double sqrtFiveOverSix = 0.9128709291752769;
    const double gamman = -2.0 * sqrtFiveOverSix * beta * std::sqrt(Sn * meff);
    const double fn_damping = -gamman * v_rel_n;
    double fn = fn_elastic + fn_damping;
    fn = std::max(fn, 0.0);

    const double G = E / (2.0 * (1.0 + nu));
    const double Geff_mindlin = G / (2.0 * (2.0 - nu));
    const double kt_ = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
    const double St = 8.0 * Geff_mindlin * std::sqrt(R_eff * overlap);
    const double gammat = -2.0 * sqrtFiveOverSix * beta * std::sqrt(St * meff);

    if(v_rel_t_mag > 1.0e-10)
    {
        history.s_t += v_rel_t * dt_contact;
        history.s_t -= normal * normal.dot(history.s_t);

        Eigen::Vector3d ft_vector = kt_ * history.s_t + gammat * v_rel_t;
        double ft_mag = ft_vector.norm();
        double ft_max = mu * fn;

        if(ft_mag > ft_max)
        {
            history.s_t *= (ft_max / ft_mag);
            ft_vector *= (ft_max / ft_mag);
        }

        force = fn * normal - ft_vector;
    }
    else
    {
        force = fn * normal;
    }

    torque = r1.cross(force);

    if(finalize || p->R52==0)
    history.t_last = p->simtime;
}
