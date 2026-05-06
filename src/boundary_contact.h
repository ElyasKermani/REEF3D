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

#ifndef BOUNDARY_CONTACT_H_
#define BOUNDARY_CONTACT_H_

#include<vector>

class lexer;
class ghostcell;
class sixdof_obj;

using namespace std;

// Base class for domain boundary contacts (ground, walls): adds loads into Xext/Yext/Zext
class boundary_contact
{
public:

    virtual ~boundary_contact() = default;

    virtual void apply(lexer *p, ghostcell *pgc, vector<sixdof_obj*> &fb_obj) = 0;
};

#endif
