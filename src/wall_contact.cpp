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

#include"wall_contact.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

#include<algorithm>
#include<cmath>
#include<iostream>

wall_contact::wall_contact()
{
    Kw_  = 3.0e6;
    mu_  = 0.3;
    xi_  = 1.0;
    nu_  = 0.01;
}

wall_contact::~wall_contact()
{
}

void wall_contact::apply(lexer *p, ghostcell*, vector<sixdof_obj*> &fb_obj)
{
    const int obj_count = static_cast<int>(fb_obj.size());

    for(int i = 0; i < obj_count; i++)
    {
        sixdof_obj* obj = fb_obj[i];

        const double x_center = obj->c_(0);
        const double y_center = obj->c_(1);
        const double radius   = obj->radius;

        const double x_left = x_center - radius;
        if(x_left < p->xcoormin)
        {
            const double overlap = p->xcoormin - x_left;
            const double k_eff   = Kw_ * radius;
            const double m_eff   = obj->Mass_fb;
            const double c_eff   = 2.0 * xi_ * std::sqrt(k_eff * m_eff);
            const double v_normal = obj->u_fb(0);
            const double F_normal = k_eff * overlap - c_eff * std::min(v_normal, 0.0);

            obj->Xext += F_normal;

            const double v_y = obj->u_fb(1);
            const double v_z = obj->u_fb(2);
            const double v_tangential = std::sqrt(v_y*v_y + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                const double F_friction_mag = mu_ * F_normal;
                const double vy_norm = v_y / std::max(nu_, v_tangential);
                const double vz_norm = v_z / std::max(nu_, v_tangential);

                obj->Yext -= F_friction_mag * std::sin(PI/2.0 * vy_norm);
                obj->Zext -= F_friction_mag * std::sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (X-min) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        const double x_right = x_center + radius;
        if(x_right > p->xcoormax)
        {
            const double overlap = x_right - p->xcoormax;

            const double k_eff   = Kw_ * radius;
            const double m_eff   = obj->Mass_fb;
            const double c_eff   = 2.0 * xi_ * std::sqrt(k_eff * m_eff);
            // Penetration increases with v_x>0; dashpot adds to the push-back (cf. x-min: min(v,0))
            const double v_normal = obj->u_fb(0);
            const double F_normal = k_eff * overlap + c_eff * std::max(v_normal, 0.0);

            obj->Xext -= F_normal;

            const double v_y = obj->u_fb(1);
            const double v_z = obj->u_fb(2);
            const double v_tangential = std::sqrt(v_y*v_y + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                const double F_friction_mag = mu_ * F_normal;
                const double vy_norm = v_y / std::max(nu_, v_tangential);
                const double vz_norm = v_z / std::max(nu_, v_tangential);

                obj->Yext -= F_friction_mag * std::sin(PI/2.0 * vy_norm);
                obj->Zext -= F_friction_mag * std::sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (X-max) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        const double y_front = y_center - radius;
        if(y_front < p->ycoormin)
        {
            const double overlap = p->ycoormin - y_front;

            const double k_eff   = Kw_ * radius;
            const double m_eff   = obj->Mass_fb;
            const double c_eff   = 2.0 * xi_ * std::sqrt(k_eff * m_eff);
            const double v_normal = obj->u_fb(1);
            const double F_normal = k_eff * overlap - c_eff * std::min(v_normal, 0.0);

            obj->Yext += F_normal;

            const double v_x = obj->u_fb(0);
            const double v_z = obj->u_fb(2);
            const double v_tangential = std::sqrt(v_x*v_x + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                const double F_friction_mag = mu_ * F_normal;
                const double vx_norm = v_x / std::max(nu_, v_tangential);
                const double vz_norm = v_z / std::max(nu_, v_tangential);

                obj->Xext -= F_friction_mag * std::sin(PI/2.0 * vx_norm);
                obj->Zext -= F_friction_mag * std::sin(PI/2.0 * vz_norm);
            }

            if(p->mpirank == 0 && p->count % 100 == 0)
            {
                cout << "Wall Contact (Y-min) - Object " << i
                     << ": overlap=" << overlap*1000.0 << " mm"
                     << " F_normal=" << F_normal << " N"
                     << endl;
            }
        }

        const double y_back = y_center + radius;
        if(y_back > p->ycoormax)
        {
            const double overlap = y_back - p->ycoormax;

            const double k_eff   = Kw_ * radius;
            const double m_eff   = obj->Mass_fb;
            const double c_eff   = 2.0 * xi_ * std::sqrt(k_eff * m_eff);
            // Penetration increases with v_y>0; dashpot adds to push-back (cf. y-min wall)
            const double v_normal = obj->u_fb(1);
            const double F_normal = k_eff * overlap + c_eff * std::max(v_normal, 0.0);

            obj->Yext -= F_normal;

            const double v_x = obj->u_fb(0);
            const double v_z = obj->u_fb(2);
            const double v_tangential = std::sqrt(v_x*v_x + v_z*v_z);

            if(v_tangential > 1.0e-10)
            {
                const double F_friction_mag = mu_ * F_normal;
                const double vx_norm = v_x / std::max(nu_, v_tangential);
                const double vz_norm = v_z / std::max(nu_, v_tangential);

                obj->Xext -= F_friction_mag * std::sin(PI/2.0 * vx_norm);
                obj->Zext -= F_friction_mag * std::sin(PI/2.0 * vz_norm);
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
