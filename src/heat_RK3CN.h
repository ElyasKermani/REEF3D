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

#ifndef HEAT_RK3CN_H_
#define HEAT_RK3CN_H_

#include"heat_base.h"

using namespace std;

// low-storage 3rd-order Runge-Kutta (Wray) with Crank-Nicolson diffusion in each substep (Rai & Moin)
class heat_RK3CN final : public heat_base
{
public:
    heat_RK3CN(lexer *, fdm*, ghostcell*,heat*&);
	virtual ~heat_RK3CN();

	void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*) override final;
	void ttimesave(lexer*, fdm*) override final;

private:
    void substep(lexer*, fdm*, ghostcell*, convection*, solver*, field&, field&, int);

    field4 ark1,ark2,E,Eold,Dk;
};

#endif
