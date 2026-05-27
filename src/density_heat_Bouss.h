/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef DENSITY_HEAT_BOUSS_H_
#define DENSITY_HEAT_BOUSS_H_

#include"density.h"
#include"increment.h"

class fdm;
class lexer;
class heat;

using namespace std;

class density_heat_Bouss final : public density, virtual public increment
{

public:
    density_heat_Bouss(lexer*,heat*&);
	virtual ~density_heat_Bouss();

	double roface(lexer*,fdm*,int,int,int) override final;

private:
	double H,roval,phival;
	double visc_1,visc_2,ro_1,ro_2;
	double T0_1,T0_2;
	double epsi;

	heat *pheat;
};

#endif
