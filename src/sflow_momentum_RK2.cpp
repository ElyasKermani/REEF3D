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

#include"sflow_momentum_RK2.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_convection.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sflow_fsf.h"
#include"sflow_forcing.h"
#include"ioflow.h"
#include"solver2D.h"
#include"sflow_rough_manning.h"
#include"sflow_rough_void.h"
#include"sflow_rheology_f.h"
#include"sflow_rheology_v.h"
#include"6DOF.h"

sflow_momentum_RK2::sflow_momentum_RK2(lexer *p, fdm2D *b, sflow_convection *pconvection, sflow_diffusion *ppdiff, sflow_pressure* ppressure,
                                                    solver2D *psolver, solver2D *ppoissonsolver, ioflow *pioflow, sflow_fsf *pfreesurf,sflow_forcing *ppsfdf,
                                                    sixdof *pp6dof)
                                                    :Prk1(p),Qrk1(p),wrk1(p),etark1(p)
{
	gcval_u=10;
	gcval_v=11;
    gcval_w=12;

    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;

	pconvec=pconvection;
	pdiff=ppdiff;
	ppress=ppressure;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
	pfsf=pfreesurf;
    p6dof=pp6dof;
    psfdf=ppsfdf;
    
    if(p->A218==0)
    prough = new sflow_rough_void(p);
    
    if(p->A218==1)
    prough = new sflow_rough_manning(p);
    
    
    if(p->W90==0)
    prheo = new sflow_rheology_v(p);
    
    if(p->W90==1)
    prheo = new sflow_rheology_f(p);
}

sflow_momentum_RK2::~sflow_momentum_RK2()
{
}

void sflow_momentum_RK2::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	        
    pflow->discharge2D(p,b,pgc);
    pflow->inflow2D(p,b,pgc,b->P,b->Q,b->bed,b->eta);
    pflow->inflow2D(p,b,pgc,Prk1,Qrk1,b->bed,b->eta);
	
//Step 1
//--------------------------------------------------------	
    // fsf
    pfsf->disc(p,b,pgc,b->P,b->Q,b->ws,b->eta);
    
    SLICELOOP4
    etark1(i,j) =      b->eta(i,j) 
    
                -      p->dt*(b->P(i,j)*b->hx(i,j) - b->P(i-1,j)*b->hx(i-1,j)
                       +      b->Q(i,j)*b->hy(i,j) - b->Q(i,j-1)*b->hy(i,j-1))/p->DXM;
                             
    pgc->gcsl_start4(p,etark1,gcval_eta);
    pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,etark1);
    pgc->gcsl_start4(p,etark1,gcval_eta);
    pfsf->disc(p,b,pgc,b->P,b->Q,b->ws,etark1);
    pfsf->breaking(p,b,pgc,etark1,b->eta,1.0);
    pflow->waterlevel2D(p,b,pgc,etark1);
    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,etark1,gcval_eta);
    

    // U
	starttime=pgc->timer();
	pflow->isource2D(p,b,pgc); 
	ppress->upgrad(p,b,etark1,b->eta);
	irhs(p,b,pgc,b->P,1.0);
    prough->u_source(p,b,b->P);
    prheo->u_source(p,b,b->P,b->Q);
    p6dof->isource2D(p,b,pgc);
	pconvec->start(p,b,b->P,1,b->P,b->Q);
	pdiff->diff_u(p,b,pgc,psolv,b->P,b->Q,1.0);

	SLICELOOP1
	Prk1(i,j) = b->P(i,j)
				+ p->dt*b->F(i,j);
	
	pgc->gcsl_start1(p,Prk1,gcval_u);

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
	ppress->vpgrad(p,b,etark1,b->eta);
	jrhs(p,b,pgc,b->Q,1.0);
    prough->v_source(p,b,b->Q);
    prheo->v_source(p,b,b->P,b->Q);
    p6dof->jsource2D(p,b,pgc);
	pconvec->start(p,b,b->Q,2,b->P,b->Q);
	pdiff->diff_v(p,b,pgc,psolv,b->P,b->Q,1.0);

	SLICELOOP2
	Qrk1(i,j) = b->Q(i,j)
			  + p->dt*b->G(i,j);
	
	pgc->gcsl_start2(p,Qrk1,gcval_v);
	
    p->vtime=pgc->timer()-starttime;
	
    // W
    ppress->wpgrad(p,b,etark1,b->eta);
    if(p->A214==1)
    pconvec->start(p,b,b->ws,4,b->P,b->Q);
    pdiff->diff_w(p,b,pgc,psolv,b->P,b->Q,b->ws,1.0);
    
    SLICELOOP4
	wrk1(i,j) = b->ws(i,j)
			  + p->dt*b->L(i,j);
              
    pgc->gcsl_start4(p,wrk1,12);
    
    psfdf->forcing(p,b,pgc,p6dof,0,1.0,Prk1,Qrk1,wrk1,etark1,b->hp,0);

	// press
    ppress->start(p,b,pgc,ppoissonsolv,pflow, Prk1, Qrk1, b->P, b->Q, wrk1, etark1, 1.0);

	pflow->pm_relax(p,pgc,b->press);

	pflow->um_relax(p,pgc,Prk1,b->bed,b->eta);
	pflow->vm_relax(p,pgc,Qrk1,b->bed,b->eta);
    pflow->wm_relax(p,pgc,wrk1,b->bed,b->eta);

	pgc->gcsl_start1(p,Prk1,gcval_u);
	pgc->gcsl_start2(p,Qrk1,gcval_v);
    pgc->gcsl_start4(p,wrk1,gcval_w);
        
//Step 2
//--------------------------------------------------------
    
    //fsf
    pfsf->disc(p,b,pgc,Prk1,Qrk1,wrk1,etark1);
    
    SLICELOOP4
    b->eta(i,j) = 0.5*b->eta(i,j) + 0.5*etark1(i,j)
    
                - 0.5*p->dt*(Prk1(i,j)*b->hx(i,j) - Prk1(i-1,j)*b->hx(i-1,j)
                       +     Qrk1(i,j)*b->hy(i,j) - Qrk1(i,j-1)*b->hy(i,j-1))/p->DXM;
                       
    pgc->gcsl_start4(p,b->eta,gcval_eta);   
    pfsf->depth_update(p,b,pgc,Prk1,Qrk1,wrk1,b->eta);
    pfsf->disc(p,b,pgc,Prk1,Qrk1,wrk1,b->eta);
    pfsf->breaking(p,b,pgc,b->eta,etark1,0.5);
    pflow->waterlevel2D(p,b,pgc,b->eta);
    pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,gcval_eta);

    
	// U
	starttime=pgc->timer();

	pflow->isource2D(p,b,pgc);
	ppress->upgrad(p,b,b->eta,etark1);
	irhs(p,b,pgc,Prk1,0.5);
    prough->u_source(p,b,Prk1);
    prheo->u_source(p,b,Prk1,Qrk1);
    p6dof->isource2D(p,b,pgc);
	pconvec->start(p,b,Prk1,1,Prk1,Qrk1);
	pdiff->diff_u(p,b,pgc,psolv,Prk1,Qrk1,0.5);

	SLICELOOP1
	b->P(i,j) = 0.5*b->P(i,j) + 0.5*Prk1(i,j)
				+ 0.5*p->dt*b->F(i,j);
				
	pgc->gcsl_start1(p,b->P,gcval_u);
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pflow->jsource2D(p,b,pgc);
	ppress->vpgrad(p,b,b->eta,etark1);
	jrhs(p,b,pgc,Qrk1,0.5);
    prough->v_source(p,b,Qrk1);
    prheo->v_source(p,b,Prk1,Qrk1);
    p6dof->jsource2D(p,b,pgc);
	pconvec->start(p,b,Qrk1,2,Prk1,Qrk1);
	pdiff->diff_v(p,b,pgc,psolv,Prk1,Qrk1,0.5);

	SLICELOOP2
	b->Q(i,j) = 0.5*b->Q(i,j) + 0.5*Qrk1(i,j)
			  + 0.5*p->dt*b->G(i,j);
	
	pgc->gcsl_start2(p,b->Q,gcval_v);	
	
    p->vtime+=pgc->timer()-starttime;
	
    // W     
    ppress->wpgrad(p,b,b->eta,etark1);
    if(p->A214==1)
    pconvec->start(p,b,wrk1,4,Prk1,Qrk1);
    pdiff->diff_w(p,b,pgc,psolv,Prk1,Qrk1,wrk1,0.5);
    
    SLICELOOP4
	b->ws(i,j) = 0.5*b->ws(i,j) + 0.5*wrk1(i,j)
			  + 0.5*p->dt*b->L(i,j);
              
    pgc->gcsl_start4(p,b->ws,12);
    
    psfdf->forcing(p,b,pgc,p6dof,1,0.5,b->P,b->Q,b->ws,b->eta,b->hp,1);
              
	//--------------------------------------------------------
	// pressure
	ppress->start(p,b,pgc,ppoissonsolv,pflow, b->P, b->Q, Prk1, Qrk1, b->ws, b->eta, 0.5);
	
    pflow->pm_relax(p,pgc,b->press);
    
	pflow->um_relax(p,pgc,b->P,b->bed,b->eta);
	pflow->vm_relax(p,pgc,b->Q,b->bed,b->eta);
    pflow->wm_relax(p,pgc,b->ws,b->bed,b->eta);

	pgc->gcsl_start1(p,b->P,gcval_u);
	pgc->gcsl_start2(p,b->Q,gcval_v);
    pgc->gcsl_start4(p,b->ws,gcval_w);

    pfsf->breaking_persist(p,b,pgc,b->eta,b->eta_n,1.0);
    
    SLICELOOP4
    b->eta_n(i,j) = b->eta(i,j);

    pgc->gcsl_start4(p,b->eta_n,gcval_eta);
}

void sflow_momentum_RK2::irhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
	n=0;
	if(p->D20<4)
	SLICELOOP1
	{
	b->maxF=MAX(fabs(b->F(i,j)),b->maxF);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

void sflow_momentum_RK2::jrhs(lexer *p, fdm2D *b, ghostcell *pgc, slice &f, double alpha)
{
	n=0;
	if(p->D20<4)
	SLICELOOP2
	{
	b->maxG=MAX(fabs(b->G(i,j)),b->maxG);
	b->rhsvec.V[n]=0.0;
	++n;
	}
}

