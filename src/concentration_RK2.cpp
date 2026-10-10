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

#include"concentration_RK2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_concentration.h"

concentration_RK2::concentration_RK2(lexer* p, fdm* a, ghostcell *pgc) : bc_concentration(p), concentration_io(p,a)
{
	gcval_concentration=80;
}

concentration_RK2::~concentration_RK2()
{
}

// The explicit terms of a stage are combined first and the implicit diffusion is solved afterwards with
// the stage weight (as momentum_rk). Before, C was overwritten by its diffused value in the first stage and
// reused in the final combination: the effective diffusivity was 1.25x (RK2) and 1.48x (RK3) too large.
void concentration_RK2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    field4 ark1(p);

    starttime=pgc->timer();

    auto rhs_of = [&](field &f)
    {
        clearrhs(p,a,pgc);
        pconvec->start(p,a,f,4,a->u,a->v,a->w);

        if(p->D20==1)
        pdiff->diff_scalar(p,a,pgc,psolv,f,a->visc,a->eddyv,1.0,1.0);
    };

    auto implicit = [&](field &f, double w)
    {
        if(p->D20>=2)
        {
        clearrhs(p,a,pgc);
        pdiff->diff_scalar(p,a,pgc,psolv,f,a->visc,a->eddyv,1.0,w);
        }
    };

// Step 1
    rhs_of(C);

	LOOP
	ark1(i,j,k) = C(i,j,k) + p->dt*a->L(i,j,k);

    implicit(ark1,1.0);
    bc_concentration_start(p,a,pgc,ark1);
	pgc->start4(p,ark1,gcval_concentration);

// Step 2
    rhs_of(ark1);

	LOOP
	C(i,j,k) = 0.5*C(i,j,k) + 0.5*ark1(i,j,k) + 0.5*p->dt*a->L(i,j,k);

    implicit(C,0.5);
    bc_concentration_start(p,a,pgc,C);
	pgc->start4(p,C,gcval_concentration);

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->concentrationtime=pgc->timer()-starttime;
}

void concentration_RK2::ttimesave(lexer *p, fdm* a)
{
}

void concentration_RK2::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
