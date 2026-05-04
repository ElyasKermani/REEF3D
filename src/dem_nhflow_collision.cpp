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

#include"dem_nhflow_collision.h"
#include"6DOF_collision.h"
#include"lexer.h"

dem_nhflow_collision::dem_nhflow_collision(lexer *p, ghostcell *pgc)
{
    p_collision = new sixdof_collision(p,pgc);
}

dem_nhflow_collision::~dem_nhflow_collision()
{
    delete p_collision;
}

void dem_nhflow_collision::ini(lexer*, ghostcell*)
{
}

void dem_nhflow_collision::initialize(lexer*, fdm_nhf*, ghostcell*)
{
}

void dem_nhflow_collision::start(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<sixdof_obj*> &fb_obj, int iter, bool finalize)
{
    if(fb_obj.size()<2)
    return;

    if(p->X49==1)
    start_sphere_sphere(p,d,pgc,fb_obj,iter,finalize);

    if(p->X49==2)
    start_triangle_sat(p,d,pgc,fb_obj,iter,finalize);

    if(p->X49==0 || p->X49==3)
    start_adaptive(p,d,pgc,fb_obj,iter,finalize);
}

void dem_nhflow_collision::start_sphere_sphere(lexer *p, fdm_nhf*, ghostcell *pgc, vector<sixdof_obj*> &fb_obj, int, bool)
{
    p_collision->set_detection_mode(CollisionDetectionMode::SphereOnly);

    policy.apply(p,fb_obj,*p_collision);
    p_collision->calculate_collision_forces(p,pgc,fb_obj);
}

void dem_nhflow_collision::start_triangle_sat(lexer *p, fdm_nhf*, ghostcell *pgc, vector<sixdof_obj*> &fb_obj, int, bool)
{
    p_collision->set_detection_mode(CollisionDetectionMode::TriangleSATOnly);

    policy.apply(p,fb_obj,*p_collision);
    p_collision->calculate_collision_forces(p,pgc,fb_obj);
}

void dem_nhflow_collision::start_adaptive(lexer *p, fdm_nhf*, ghostcell *pgc, vector<sixdof_obj*> &fb_obj, int, bool)
{
    p_collision->set_detection_mode(CollisionDetectionMode::Adaptive);

    policy.apply(p,fb_obj,*p_collision);
    p_collision->calculate_collision_forces(p,pgc,fb_obj);
}
