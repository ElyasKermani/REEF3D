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

#ifndef DEM_NHFLOW_COLLISION_H_
#define DEM_NHFLOW_COLLISION_H_

#include"dem_nhflow.h"
#include"dem_policy.h"

class lexer;
class fdm_nhf;
class ghostcell;
class sixdof_obj;
class dem_collision;

using namespace std;

class dem_nhflow_collision : public dem_nhflow
{
public:

    dem_nhflow_collision(lexer*, ghostcell*);
    virtual ~dem_nhflow_collision();

    void ini(lexer*, ghostcell*) override;
    void initialize(lexer*, fdm_nhf*, ghostcell*) override;
    void start(lexer*, fdm_nhf*, ghostcell*, vector<sixdof_obj*>&, int, bool) override;

    void start_sphere_sphere(lexer*, fdm_nhf*, ghostcell*, vector<sixdof_obj*>&, int, bool);
    void start_triangle_sat(lexer*, fdm_nhf*, ghostcell*, vector<sixdof_obj*>&, int, bool);
    void start_adaptive(lexer*, fdm_nhf*, ghostcell*, vector<sixdof_obj*>&, int, bool);

private:

    dem_collision *p_collision;
    dem_policy policy;
};

#endif
