/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/
#include"ptf_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"ptf_laplace_cds2.h"
#include"ptf_laplace_cds4.h"
#include"onephase.h"
#include"ptf_fsf_update.h"
#include"ptf_bed_update.h"

ptf_RK3::ptf_RK3(lexer *p, fdm *a, ghostcell *pgc) : ptf_fsfbc(p,a,pgc),erk1(p),erk2(p),frk1(p),frk2(p)
{
    gcval=250;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    gcval_eta = 50;
    gcval_fifsf = 50;
    
    if(p->A320==1)
    plap = new ptf_laplace_cds2;
    
    if(p->A320==2)
    plap = new ptf_laplace_cds4;
    
    pfsfupdate = new ptf_fsf_update(p,a,pgc);
    
    pbedupdate = new ptf_bed_update(p,a,pgc);
}

ptf_RK3::~ptf_RK3()
{
}

void ptf_RK3::start(lexer *p, fdm *a, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	
// Step 1
    fsfdisc(p,a,pgc,a->eta,a->Fifsf,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);

    SLICELOOP4
	erk1(i,j) = a->eta(i,j) + p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,erk1,gcval_eta);

    // fsf Fi
    dfsfbc(p,a,pgc,a->eta);

    SLICELOOP4
	frk1(i,j) = a->Fifsf(i,j) + p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,frk1,gcval_fifsf);

    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk1);
    pflow->fifsf_relax(p,pgc,frk1);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,erk1);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk1,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);

// Step 2
    fsfdisc(p,a,pgc,erk1,frk1,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
	erk2(i,j) = 0.75*a->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,erk2,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk1);
    
    SLICELOOP4
	frk2(i,j) = 0.75*a->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,frk2,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,erk2);
    pflow->fifsf_relax(p,pgc,frk2);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,erk2);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,frk2,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);

// Step 3
    fsfdisc(p,a,pgc,erk2,frk2,a->Fi);
    
    // fsf eta
    kfsfbc(p,a,pgc);
    
    SLICELOOP4
	a->eta(i,j) = (1.0/3.0)*a->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,a->eta,gcval_eta);
    
    // fsf Fi
    dfsfbc(p,a,pgc,erk2);
    
    SLICELOOP4
	a->Fifsf(i,j) = (1.0/3.0)*a->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*a->K(i,j);
    
    pgc->gcsl_start4(p,a->Fifsf,gcval_fifsf);
    
    // Set Boundary Conditions
    pflow->eta_relax(p,pgc,a->eta);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    pbedupdate->waterdepth(p,a,pgc);
    pbedupdate->bedbc(p,a,pgc,a->Fi);
    
    // solve Fi
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,gcval);
    plap->start(p,a,pgc,psolv,a->Fi);

    FLUIDLOOP
    a->test(i,j,k) = a->Fifsf(i,j);

    pfsfupdate->velcalc(p,a,pgc,a->Fi);
}

void ptf_RK3::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reini *preini, onephase *poneph)
{	
    pfsfupdate->fsfupdate(p,a,pgc,pflow,poneph,a->eta);
    pfsfupdate->etaloc(p,a,pgc);
    
    // potential ini
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pflow->fifsf_relax(p,pgc,a->Fifsf);
    pgc->start4(p,a->Fi,250);
    pgc->gcsl_start4(p,a->Fifsf,50);
    
    pfsfupdate->fsfbc(p,a,pgc,a->Fifsf,a->Fi);
    
    pfsfupdate->fsfepol(p,a,pgc,a->eta,a->Fi);
    
    
    pgc->gcsl_start4(p,a->eta,50);
}

void ptf_RK3::inidisc(lexer *p, fdm *a, ghostcell *pgc)
{	
}
