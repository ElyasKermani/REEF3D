/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef BC_CONCENTRATION_H_
#define BC_CONCENTRATION_H_

#include"increment.h"

class lexer;
class fdm;
class field;
class ghostcell;

using namespace std;

class bc_concentration : public increment
{
public:
	bc_concentration(lexer*);
	virtual ~bc_concentration();
	void bc_concentration_start(lexer*,fdm*,ghostcell*,field&);

private:


};
#endif


