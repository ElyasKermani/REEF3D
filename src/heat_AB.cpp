/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
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
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_heat.h"

heat_AB::heat_AB(lexer* p, fdm* a, ghostcell *pgc, heat *&pheat) : bcheat(p), heat_print(p,a), tab(p), thermdiff(p)
{
	gcval_heat=80;
}

heat_AB::~heat_AB()
{
}

void heat_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{

    starttime=pgc->timer();
    
    diff_update(p,a,pgc);

    clearrhs(p,a,pgc);
	pconvec->start(p,a,T,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,T,thermdiff,a->eddyv,p->sigT,1.0);

	if(p->count==1)
	LOOP
	tab(i,j,k)=a->L(i,j,k);

	LOOP
	{
	T(i,j,k)+=p->dt*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->L(i,j,k)
								-(p->dt/p->dt_old)*tab(i,j,k));
	tab(i,j,k)=a->L(i,j,k);
	}
	
    bcheat_start(p,a,pgc,T);
	pgc->start4(p,T,gcval_heat);
	pupdate->start(p,a,pgc);

	p->heattime=pgc->timer()-starttime;
}

void heat_AB::ttimesave(lexer *p, fdm* a)
{
}

void heat_AB::diff_update(lexer *p, fdm *a, ghostcell *pgc)
{
    double alpha_1;
	double alpha_2;
    double H;
    double epsi=p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
    
    if(p->H9==1)
    {
    alpha_1 = p->H1;
	alpha_2 = p->H2;
    }
    
    if(p->H9==2)
    {
    alpha_1 = p->H2;
	alpha_2 = p->H1;
    }
    
    
    
    LOOP
	{
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		thermdiff(i,j,k) = alpha_1*H + alpha_2*(1.0-H);
	}
    
    pgc->start4(p,thermdiff,1);
}

void heat_AB::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
