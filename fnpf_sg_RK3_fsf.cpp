/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"fnpf_sg_RK3_fsf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_sg_laplace_cds2.h"
#include"fnpf_sg_laplace_cds2_v2.h"
#include"fnpf_sg_laplace_cds4.h"
#include"fnpf_sg_laplace_cds4_bc2.h"
#include"fnpf_sg_laplace_cds24.h"
#include"fnpf_sg_laplace_HOS.h"
#include"onephase.h"
#include"fnpf_sg_fsfbc.h"
#include"fnpf_sg_fsfbc_wd.h"

fnpf_sg_RK3_fsf::fnpf_sg_RK3_fsf(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_sg_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      erk1(p),erk2(p),frk1(p),frk2(p)
{
    gcval=250;
    if(p->j_dir==0)
    gcval=150;
   
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
    
    
    if(p->A320==1)
    plap = new fnpf_sg_laplace_cds2(p);
    
    if(p->A320==2)
    plap = new fnpf_sg_laplace_cds4(p);
    
    if(p->A320==3)
    plap = new fnpf_sg_laplace_cds4_bc2(p);
    
    if(p->A320==4)
    plap = new fnpf_sg_laplace_cds24(p);
    
    if(p->A320==5)
    plap = new fnpf_sg_laplace_cds2_v2(p);
    
    
    if(p->A320==11)
    plap = new fnpf_sg_laplace_HOS(p);
    
    
    if(p->A343==0)
    pf = new fnpf_sg_fsfbc(p,c,pgc);
    
    if(p->A343>=1)
    pf = new fnpf_sg_fsfbc_wd(p,c,pgc);
    
    t0=0.0;
}

fnpf_sg_RK3_fsf::~fnpf_sg_RK3_fsf()
{
}

void fnpf_sg_RK3_fsf::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini, onephase* poneph)
{	   
    
    LOOP
    c->test(i,j,k)=0.0;
    
// Step 1
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    integrate(p,c,pgc);
    
    SLICELOOP4
    erk1(i,j) = c->eta(i,j) - p->dt*((c->P(i,j) - c->P(i-1,j))/p->DXN[IP] + (c->Q(i,j)- c->Q(i,j-1))/p->DYN[JP]);
    
    pf->damping(p,c,pgc,erk1,gcval_eta,1.0);
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,c->eta);

    SLICELOOP4
	frk1(i,j) = c->Fifsf(i,j) + p->dt*c->K(i,j);
    
    pf->damping(p,c,pgc,frk1,gcval_fifsf,1.0);
   
    pflow->eta_relax(p,pgc,erk1);
    pgc->gcsl_start4(p,erk1,gcval_eta);
    pf->coastline(p,c,pgc,erk1);
    pf->coastline(p,c,pgc,frk1);
    pflow->fifsf_relax(p,pgc,frk1);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk1, c->eta, frk1,1.0);
    pf->wetdry(p,c,pgc,erk1,frk1);
    pf->fsfdisc(p,c,pgc,erk1,frk1);
    sigma_update(p,c,pgc,pf,erk1);
  
    // Set Boundary Conditions Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk1,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,erk1,frk1);
    
     LOOP
    if(c->breaking(i,j)==1)
    c->test(i,j,k)=1.0;

// Step 2
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,frk1,erk1);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    integrate(p,c,pgc);
    
    SLICELOOP4
    erk2(i,j) = 0.75*c->eta(i,j)  + 0.25*erk1(i,j) - 0.25*p->dt*((c->P(i,j) - c->P(i-1,j))/p->DXN[IP] + (c->Q(i,j)-c->Q(i,j-1))/p->DYN[JP]);
    
    pf->damping(p,c,pgc,erk2,gcval_eta,0.25);
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,erk1);
    
    SLICELOOP4
	frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*c->K(i,j);
    
    pf->damping(p,c,pgc,frk2,gcval_fifsf,0.25);
    
    pflow->eta_relax(p,pgc,erk2);
    pgc->gcsl_start4(p,erk2,gcval_eta);
    pf->coastline(p,c,pgc,erk2);
    pf->coastline(p,c,pgc,frk2);
    pflow->fifsf_relax(p,pgc,frk2);
    pgc->gcsl_start4(p,frk2,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, erk2, erk1, frk2, 0.25);
    pf->wetdry(p,c,pgc,erk2,frk2);
    pf->fsfdisc(p,c,pgc,erk2,frk2);
    sigma_update(p,c,pgc,pf,erk2);
    
    // Set Boundary Conditions Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,frk2,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,erk2,frk2);
    
     LOOP
    if(c->breaking(i,j)==1)
    c->test(i,j,k)=1.0;

// Step 3 
    pflow->inflow_fnpf(p,pgc,c->Fi,c->Uin,frk2,erk2);
    
    // fsf eta
    pf->kfsfbc(p,c,pgc);
    
    integrate(p,c,pgc);
    
    SLICELOOP4
    c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) - (2.0/3.0)*p->dt*((c->P(i,j) - c->P(i-1,j))/p->DXN[IP] + (c->Q(i,j)-c->Q(i,j-1))/p->DYN[JP]);
    
    pf->damping(p,c,pgc,c->eta,gcval_eta,2.0/3.0);
    
    // fsf Fi
    pf->dfsfbc(p,c,pgc,erk2);
    
    SLICELOOP4
	c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*c->K(i,j);
    
    pf->damping(p,c,pgc,c->Fifsf,gcval_fifsf,2.0/3.0);
    
    pflow->eta_relax(p,pgc,c->eta);
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pf->coastline(p,c,pgc,c->eta);
    pf->coastline(p,c,pgc,c->Fifsf);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    // fsfdisc and sigma update
    pf->breaking(p, c, pgc, c->eta, erk2,c->Fifsf,2.0/3.0);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    // Set Boundary Conditions Fi
    pflow->fivec_relax(p,pgc,c->Fi);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    // solve Fi
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap->start(p,c,pgc,psolv,pf,c->Fi);
    pflow->fivec_relax(p,pgc,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    

    //---------------------------------
    //LOOP
    //c->test(i,j,k)=c->vb(i,j);
    
    LOOP
    if(c->breaking(i,j)==1)
    c->test(i,j,k)=1.0;
    
    pgc->start4(p,c->test,50);

    bedbc_sig(p,c,pgc,c->Fi,pf);
    velcalc_sig(p,c,pgc,c->Fi);
    
    t0 = p->simtime;
}

void fnpf_sg_RK3_fsf::inidisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{	
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    etaloc_sig(p,c,pgc);
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    sigma_ini(p,c,pgc,pf,c->eta);
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);

    pgc->start4(p,c->test,50);
    
    
    pf->coastline(p,c,pgc,c->eta);
    pf->coastline(p,c,pgc,c->Fifsf);
    
    
    velcalc_sig(p,c,pgc,c->Fi);
    
    pgc->start7V(p,c->U,c->bc,210);
    pgc->start7V(p,c->V,c->bc,210);
    pgc->start7V(p,c->W,c->bc,210);
    
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
}

void fnpf_sg_RK3_fsf::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini

    pf->coastline(p,c,pgc,c->eta);
    pf->coastline(p,c,pgc,c->Fifsf);
}

void fnpf_sg_RK3_fsf::integrate(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    double uvel,vvel,sigz;
    
    SLICELOOP1
    c->P(i,j) = 0.0;
    
    SLICELOOP2
    c->Q(i,j) = 0.0;
    
    SLICELOOP1
    FKLOOP
    FPCHECK
    if(k>0)
    {
        uvel = (c->Fi[FIp1JK]-c->Fi[FIJK])/p->DXP[IP];
        
        if(k<p->knoz)
        uvel += 0.5*(p->sigx[FIJK]+p->sigx[FIp1JK])*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
        
        if(k==p->knoz)
        uvel += 0.5*(p->sigx[FIJK]+p->sigx[FIp1JK])*((c->Fi[FIJK]-c->Fi[FIJKm1])/p->DZN[KP]);

        c->P(i,j) += uvel*p->DZP[KP]*0.5*(p->sigz[IJ]+p->sigz[Ip1J]);
    }
    
    if(p->j_dir==1)
    SLICELOOP2
    FKLOOP
    FPCHECK
    if(k>0)
    {
        vvel = (c->Fi[FIJp1K]-c->Fi[FIJK])/p->DYP[JP];
        
        if(k<p->knoz)
        vvel += 0.5*(p->sigy[FIJK]+p->sigy[FIJp1K])*((c->Fi[FIJKp1]-c->Fi[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));
        
        if(k==p->knoz)
        vvel += 0.5*(p->sigy[FIJK]+p->sigy[FIJp1K])*((c->Fi[FIJK]-c->Fi[FIJKm1])/p->DZN[KP]);

        c->Q(i,j) += vvel*p->DZN[KP]*0.5*(p->sigz[IJ]+p->sigz[IJp1]);
    }
    
    
	pgc->gcsl_start1(p,c->P,10);
    pgc->gcsl_start2(p,c->Q,11);
}






