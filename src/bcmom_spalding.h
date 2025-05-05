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

#ifndef BCMOM_SPALDING_H_
#define BCMOM_SPALDING_H_

#include"surftens.h"
#include"wall_functions.h"

class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;

using namespace std;

class bcmom_spalding : public surftens, public wall_function_spalding
{
public:
	bcmom_spalding(lexer*);
	virtual ~bcmom_spalding();
	virtual void bcmom_start(fdm*, lexer*, ghostcell*, turbulence*, field&, int);

private:
	const double kappa;
	int gcval_phi, bckin;
};

#endif 