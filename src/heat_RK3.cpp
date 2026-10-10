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

#include"heat_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"fluid_update.h"

heat_RK3::heat_RK3(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : heat_base(p,a), ark1(p), ark2(p), E(p), D1(p), D2(p)
{
}

heat_RK3::~heat_RK3()
{
}

/*--------------------------------------------------------------------
explicit part: Shu-Osher TVD RK3
implicit diffusion (dmode 2), with D the same discrete operator:
    (I -     dt D) T1     = T + dt E(T)
    (I - 1/4 dt D) T2     = 3/4 T + 1/4 T1 + 1/4 dt E(T1)
    (I - 1/2 dt D) T(n+1) = 1/3 T + 2/3 T2 + 2/3 dt E(T2) + dt (-2/3 D(T1) + 5/6 D(T2))
implicit weights (-1/2, 1, 1/2) on (T1, T2, T(n+1)) and (1/4, 1/4) for T2:
2nd order in time together with the explicit part and L-stable
(same combination as the IMEX option of momentum_rk)
--------------------------------------------------------------------*/

void heat_RK3::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    diff_update(p,a,pgc);

// Step 1
    heat_E(p,a,pgc,pconvec,T,E);

	LOOP
	ark1(i,j,k) = T(i,j,k) + p->dt*E(i,j,k);

    if(dmode==2)
    heat_solve(p,a,pgc,psolv,ark1,1.0);

    heat_bc(p,a,pgc,ark1);

// Step 2
    heat_E(p,a,pgc,pconvec,ark1,E);

	LOOP
	ark2(i,j,k) = 0.75*T(i,j,k) + 0.25*ark1(i,j,k) + 0.25*p->dt*E(i,j,k);

    if(dmode==2)
    {
    heat_D(p,a,ark1,D1);
    heat_solve(p,a,pgc,psolv,ark2,0.25);
    }

    heat_bc(p,a,pgc,ark2);

// Step 3
    heat_E(p,a,pgc,pconvec,ark2,E);

	LOOP
	T(i,j,k) = (1.0/3.0)*T(i,j,k) + (2.0/3.0)*ark2(i,j,k) + (2.0/3.0)*p->dt*E(i,j,k);

    if(dmode==2)
    {
        heat_D(p,a,ark2,D2);

        LOOP
        T(i,j,k) += p->dt*(-(2.0/3.0)*D1(i,j,k) + (5.0/6.0)*D2(i,j,k));

        heat_solve(p,a,pgc,psolv,T,0.5);
    }

    heat_bc(p,a,pgc,T);

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->heattime=pgc->timer()-starttime;
}

void heat_RK3::ttimesave(lexer *p, fdm* a)
{
}
