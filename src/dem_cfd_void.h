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

#ifndef DEM_CFD_VOID_H_
#define DEM_CFD_VOID_H_

#include"dem_cfd.h"

class lexer;
class fdm;
class ghostcell;
class sixdof_obj;

using namespace std;

class dem_cfd_void : public dem_cfd
{
public:

    dem_cfd_void();
    virtual ~dem_cfd_void();

    void ini(lexer*, ghostcell*) override;
    void initialize(lexer*, fdm*, ghostcell*) override;
    void start(lexer*, fdm*, ghostcell*, vector<sixdof_obj*>&, int, bool, double dt_contact) override;
};

#endif
