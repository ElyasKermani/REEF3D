/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef BC_IKEPSILON_SPALDING_H_
#define BC_IKEPSILON_SPALDING_H_

#include"increment.h"
#include"wall_functions.h"

class fdm;
class lexer;
class field;
class turbulence;

using namespace std;

class bc_ikepsilon_spalding : public increment, public wall_function_spalding
{
public:
	bc_ikepsilon_spalding(lexer*);
	virtual ~bc_ikepsilon_spalding();
	void bckeps_start(fdm*, lexer*, field&, field&, int, turbulence*);

private:
	int count,q;
	double fac,value;
};

#endif 