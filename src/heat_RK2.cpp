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

#include"heat_RK2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"fluid_update.h"

heat_RK2::heat_RK2(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : heat_base(p,a), ark1(p), E(p), D0(p), D1(p)
{
}

heat_RK2::~heat_RK2()
{
}

/*--------------------------------------------------------------------
explicit part: T1 = T + dt E(T),  T(n+1) = 1/2 T + 1/2 (T1 + dt E(T1))
implicit diffusion (dmode 2), with D the same discrete operator:
    (I - dt D) T1      = T + dt E(T)
    (I - dt D) T(n+1)  = 1/2 T + 1/2 T1 + 1/2 dt E(T1) + dt (1/2 D(T) - D(T1))
implicit weights (1/2, -1/2, 1) on (T, T1, T(n+1)): 2nd order in time
together with the explicit part, A-stable, |G| -> 1/2 for stiff modes
--------------------------------------------------------------------*/

void heat_RK2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    diff_update(p,a,pgc);

// Step 1
    heat_E(p,a,pgc,pconvec,T,E);

    if(dmode==2)
    heat_D(p,a,T,D0);

	LOOP
	ark1(i,j,k) = T(i,j,k) + p->dt*E(i,j,k);

    if(dmode==2)
    heat_solve(p,a,pgc,psolv,ark1,1.0);

    heat_bc(p,a,pgc,ark1);

// Step 2
    heat_E(p,a,pgc,pconvec,ark1,E);

	LOOP
	T(i,j,k) = 0.5*T(i,j,k) + 0.5*ark1(i,j,k) + 0.5*p->dt*E(i,j,k);

    if(dmode==2)
    {
        heat_D(p,a,ark1,D1);

        LOOP
        T(i,j,k) += p->dt*(0.5*D0(i,j,k) - D1(i,j,k));

        heat_solve(p,a,pgc,psolv,T,1.0);
    }

    heat_bc(p,a,pgc,T);

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->heattime=pgc->timer()-starttime;
}

void heat_RK2::ttimesave(lexer *p, fdm* a)
{
}
