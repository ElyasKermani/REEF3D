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

#include"ground_contact.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

#include<algorithm>
#include<cmath>
#include<iostream>

ground_contact::ground_contact()
{
    zg_  = 0.0;
    Kg_  = 3.0e6;
    mu_  = 0.3;
    xi_  = 1.0;
    nu_  = 0.01;
}

ground_contact::~ground_contact()
{
}

void ground_contact::apply(lexer *p, ghostcell*, vector<sixdof_obj*> &fb_obj)
{
    const int obj_count = static_cast<int>(fb_obj.size());

    for(int i = 0; i < obj_count; i++)
    {
        sixdof_obj* obj = fb_obj[i];

        const double z_bottom = obj->c_(2) - obj->radius;

        if((zg_ - z_bottom) <= 0.0)
        {
            continue;
        }

        const double overlap = zg_ - z_bottom;
        const double k_eff   = Kg_ * obj->radius;
        const double m_eff   = obj->Mass_fb;
        const double c_eff   = 2.0 * xi_ * std::sqrt(k_eff * m_eff);
        const double v_normal = obj->u_fb(2);
        const double F_normal = k_eff * overlap - c_eff * std::max(v_normal, 0.0);

        obj->Zext_dem += F_normal;

        const double v_x = obj->u_fb(0);
        const double v_y = obj->u_fb(1);
        const double v_horiz = std::sqrt(v_x*v_x + v_y*v_y);

        const double vx_norm = v_x / std::max(nu_, v_horiz);
        const double vy_norm = v_y / std::max(nu_, v_horiz);

        const double F_friction_mag = mu_ * F_normal;
        const double Fx_friction = -F_friction_mag * std::sin(PI/2.0 * vx_norm);
        const double Fy_friction = -F_friction_mag * std::sin(PI/2.0 * vy_norm);

        obj->Xext_dem += Fx_friction;
        obj->Yext_dem += Fy_friction;

        if(p->mpirank == 0 && p->count % 100 == 0)
        {
            cout << "Ground Contact - Object " << i
                 << ": z_bottom=" << z_bottom
                 << " overlap=" << overlap*1000.0 << " mm"
                 << " F_normal=" << F_normal << " N"
                 << " F_friction=" << std::sqrt(Fx_friction*Fx_friction + Fy_friction*Fy_friction) << " N"
                 << " v_z=" << v_normal << " m/s"
                 << endl;
        }
    }
}
