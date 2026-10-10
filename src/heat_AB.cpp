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

#include"heat_AB.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"fluid_update.h"

heat_AB::heat_AB(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : heat_base(p,a), E(p), Eold(p), Dn(p)
{
    first=true;
    dt_prev=0.0;
}

heat_AB::~heat_AB()
{
}

void heat_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    diff_update(p,a,pgc);

    heat_E(p,a,pgc,pconvec,T,E);

    if(dmode==2)
    heat_D(p,a,T,Dn);

    // first step: forward Euler for the explicit part
    // (the step-size ratio uses the heat solver's own previous step: p->dt_old is not set by every time stepping)
    if(first)
    LOOP
    Eold(i,j,k)=E(i,j,k);

    const double r = first ? 1.0 : p->dt/dt_prev;

    LOOP
    {
    T(i,j,k) += p->dt*((1.0+0.5*r)*E(i,j,k) - 0.5*r*Eold(i,j,k));
    Eold(i,j,k) = E(i,j,k);
    }

    if(dmode==2)
    {
        LOOP
        T(i,j,k) += 0.5*p->dt*Dn(i,j,k);

        heat_solve(p,a,pgc,psolv,T,0.5);
    }

    heat_bc(p,a,pgc,T);

    first=false;
    dt_prev=p->dt;

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->heattime=pgc->timer()-starttime;
}

void heat_AB::ttimesave(lexer *p, fdm* a)
{
}
