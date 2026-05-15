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

#include"dem_policy.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include<algorithm>

namespace
{
    // Triangle-count thresholds for the adaptive narrow-phase.
    // Below SIMPLE_TRI_MAX  : sphere–sphere branch (fast, coarse).
    // Below MODERATE_TRI_MAX: SAT triangle–triangle without BVH.
    // Above MODERATE_TRI_MAX: SAT with body-frame BVH pruning.
    constexpr int SIMPLE_TRI_MAX   = 50;
    constexpr int MODERATE_TRI_MAX = 500;
}

dem_policy::dem_policy()
{
}

dem_config dem_policy::select(lexer *p, const vector<sixdof_obj*> &fb_obj) const
{
    dem_config config;

    // R12 now selects the detection mode at the dispatcher level (see dem_*_collision.cpp);
    // when start_adaptive is reached, adaptive routing is always desired.
    config.adaptive     = true;
    config.simple_tri   = SIMPLE_TRI_MAX;
    config.moderate_tri = MODERATE_TRI_MAX;

    // Contact-force model (R11):
    //  0 -- automatic selection (default)
    //  1 -- Linear spring-dashpot
    //  2 -- Hertz
    //  3 -- Hertz-Mindlin
    if(p->R11==1)
    {
        config.contact_model = ContactForceModel::Linear;
        return config;
    }

    if(p->R11==2)
    {
        config.contact_model = ContactForceModel::Hertz;
        return config;
    }

    if(p->R11==3)
    {
        config.contact_model = ContactForceModel::HertzMindlin;
        return config;
    }

    // Automatic selection: pick an extended DEM contact variant based on object count and size
    double max_radius = 0.0;
    for(size_t n=0; n<fb_obj.size(); ++n)
    max_radius = max(max_radius, fb_obj[n]->radius);

    if(fb_obj.size()<=2 && max_radius<0.5)
    config.contact_model = ContactForceModel::DemLinearLimited;
    else
    config.contact_model = ContactForceModel::DemHertzianNonLimited;

    return config;
}

void dem_policy::apply(lexer *p, const vector<sixdof_obj*> &fb_obj, dem_collision &collision) const
{
    const dem_config config = select(p,fb_obj);

    collision.set_contact_force_model(config.contact_model);
    collision.set_adaptive_detection(config.adaptive, config.simple_tri, config.moderate_tri);
}
