/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"kepsilon_IM1_spalding.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"convection.h"
#include"turbulence.h"

kepsilon_IM1_spalding::kepsilon_IM1_spalding(lexer* p, fdm* a, ghostcell *pgc) : 
    kepsilon_IM1(p, a, pgc),
    bc_ikep_spalding(p)
{
}

kepsilon_IM1_spalding::~kepsilon_IM1_spalding()
{
}

void kepsilon_IM1_spalding::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, vrans *pvrans)
{
    // Create a pointer to this turbulence model for use with the wall functions
    turbulence *pturb = dynamic_cast<turbulence*>(this);
    
	wallf_update(p,a,pgc,wallf);

// kin
    starttime=pgc->timer();
	clearrhs(p,a);
    pconvec->start(p,a,kin,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,kin,a->eddyv,ke_sigma_k,1.0);
	kinsource(p,a,pvrans);
	timesource(p,a,kn);
    // Use Spalding wall function for k
    bc_ikep_spalding.bckeps_start(a,p,kin,eps,gcval_kin,pturb);
	psolv->start(p,a,pgc,kin,a->rhsvec,4);
	pgc->start4(p,kin,gcval_kin);
	p->kintime=pgc->timer()-starttime;
	p->kiniter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kiniter: "<<p->kiniter<<"  kintime: "<<setprecision(3)<<p->kintime<<endl;

// eps
    starttime=pgc->timer();
	clearrhs(p,a);
    pconvec->start(p,a,eps,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,eps,a->eddyv,ke_sigma_e,1.0);
	epssource(p,a,pvrans);
	timesource(p,a,en);
	psolv->start(p,a,pgc,eps,a->rhsvec,4);
	epsfsf(p,a,pgc);
    // Use Spalding wall function for epsilon
    bc_ikep_spalding.bckeps_start(a,p,kin,eps,gcval_eps,pturb);
	pgc->start4(p,eps,gcval_eps);
	p->epstime=pgc->timer()-starttime;
	p->epsiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"epsiter: "<<p->epsiter<<"  epstime: "<<setprecision(3)<<p->epstime<<endl;

	eddyvisc(a,p,pgc,pvrans);
	pgc->start4(p,a->eddyv,24);
} 