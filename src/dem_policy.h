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

#ifndef DEM_POLICY_H_
#define DEM_POLICY_H_

#include"dem_collision.h"
#include<vector>

class lexer;
class sixdof_obj;

using namespace std;

struct dem_config
{
    ContactForceModel contact_model;  // object–object contact law
    bool adaptive;        // enable adaptive narrow-phase (triangle-count rules)
    int simple_tri;       // max triangle count for sphere–sphere branch
    int moderate_tri;     // max triangle count for SAT without BVH; above → BVH path

    dem_config()
    : contact_model(ContactForceModel::Linear),
      adaptive(true),
      simple_tri(50),
      moderate_tri(500) {}
};

class dem_policy
{
public:

    dem_policy();

    dem_config select(lexer*, const vector<sixdof_obj*>&) const;
    void apply(lexer*, const vector<sixdof_obj*>&, dem_collision&) const;
};

#endif
