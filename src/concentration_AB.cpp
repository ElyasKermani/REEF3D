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

#include"concentration_AB.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_concentration.h"

concentration_AB::concentration_AB(lexer* p, fdm* a, ghostcell *pgc) : bc_concentration(p), concentration_io(p,a), cab(p)
{
	gcval_concentration=80;
    first=true;
    dt_prev=0.0;
}

concentration_AB::~concentration_AB()
{
}

void concentration_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();

    // explicit right-hand side of a stage: convection and, with explicit diffusion (D 20 1), diffusion
    auto rhs_of = [&](field &f)
    {
        clearrhs(p,a,pgc);
        pconvec->start(p,a,f,4,a->u,a->v,a->w);

        if(p->D20==1)
        pdiff->diff_scalar(p,a,pgc,psolv,f,a->visc,a->eddyv,1.0,1.0);
    };

    // implicit diffusion (D 20 >= 2) of the combined stage with its weight: f = (I - w dt D)^-1 f
    auto implicit = [&](field &f, double w)
    {
        if(p->D20>=2)
        {
        clearrhs(p,a,pgc);
        pdiff->diff_scalar(p,a,pgc,psolv,f,a->visc,a->eddyv,1.0,w);
        }
    };

    rhs_of(C);

    // AB2 with its own previous step size (p->dt_old is not set by the fixed time stepping), Euler first
    if(first)
	LOOP
	cab(i,j,k)=a->L(i,j,k);

    const double r = first ? 1.0 : p->dt/dt_prev;

	LOOP
	{
	C(i,j,k) += p->dt*((1.0+0.5*r)*a->L(i,j,k) - 0.5*r*cab(i,j,k));
	cab(i,j,k)=a->L(i,j,k);
	}

    implicit(C,1.0);
    bc_concentration_start(p,a,pgc,C);
	pgc->start4(p,C,gcval_concentration);

    first=false;
    dt_prev=p->dt;

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->concentrationtime=pgc->timer()-starttime;
}

void concentration_AB::ttimesave(lexer *p, fdm* a)
{
}

void concentration_AB::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
