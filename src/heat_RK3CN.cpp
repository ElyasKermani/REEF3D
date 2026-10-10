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

#include"heat_RK3CN.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"solver.h"
#include"fluid_update.h"

heat_RK3CN::heat_RK3CN(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : heat_base(p,a), ark1(p), ark2(p), E(p), Eold(p), Dk(p)
{
}

heat_RK3CN::~heat_RK3CN()
{
}

/*--------------------------------------------------------------------
substep k = 0,1,2:
    (I - a_k dt D) T(k+1) = T(k) + dt (g_k E(T(k)) + z_k E(T(k-1))) + a_k dt D(T(k))
g = (8/15, 5/12, 3/4), z = (0, -17/60, -5/12), a_k = (g_k + z_k)/2
3rd order for the explicit part, 2nd order (Crank-Nicolson) for implicit diffusion.
With explicit diffusion (D 20 = 1) the diffusion is part of E and no system is solved.
--------------------------------------------------------------------*/

void heat_RK3CN::substep(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, solver *psolv, field &fin, field &fout, int s)
{
    const double g[3] = {8.0/15.0, 5.0/12.0, 3.0/4.0};
    const double z[3] = {0.0, -17.0/60.0, -5.0/12.0};
    const double ak = 0.5*(g[s]+z[s]);

    heat_E(p,a,pgc,pconvec,fin,E);

    if(s==0)
    LOOP
    Eold(i,j,k)=0.0;

    if(dmode==2)
    heat_D(p,a,fin,Dk);

    LOOP
    {
    fout(i,j,k) = fin(i,j,k) + p->dt*(g[s]*E(i,j,k) + z[s]*Eold(i,j,k));
    Eold(i,j,k) = E(i,j,k);
    }

    if(dmode==2)
    {
        LOOP
        fout(i,j,k) += ak*p->dt*Dk(i,j,k);

        heat_solve(p,a,pgc,psolv,fout,ak);
    }

    heat_bc(p,a,pgc,fout);
}

void heat_RK3CN::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    diff_update(p,a,pgc);

    substep(p,a,pgc,pconvec,psolv,T,ark1,0);
    substep(p,a,pgc,pconvec,psolv,ark1,ark2,1);
    substep(p,a,pgc,pconvec,psolv,ark2,T,2);

    pupdate->start(p,a,pgc,a->u,a->v,a->w);

    p->heattime=pgc->timer()-starttime;
}

void heat_RK3CN::ttimesave(lexer *p, fdm* a)
{
}
