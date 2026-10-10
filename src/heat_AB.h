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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef HEAT_AB_H_
#define HEAT_AB_H_

#include"heat_base.h"

using namespace std;

// 2nd-order Adams-Bashforth for convection (and explicit diffusion), Crank-Nicolson for implicit diffusion
class heat_AB final : public heat_base
{
public:
    heat_AB(lexer *, fdm*, ghostcell*,heat*&);
	virtual ~heat_AB();

	void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*) override final;
	void ttimesave(lexer*, fdm*) override final;

private:
	field4 E,Eold,Dn;
    bool first;
    double dt_prev;
};

#endif
