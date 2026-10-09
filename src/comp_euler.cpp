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

#include"comp_euler.h"
#include"lexer.h"
#include"ghostcell.h"
#include"field4.h"
#include"lagoon_output.h"
#include<cmath>
#include<iostream>
#include<algorithm>
#include<string>

using namespace std;

comp_euler::comp_euler(lexer *pp, ghostcell *ppgc) : p(pp), pgc(ppgc)
{
    for(int q=0; q<5; ++q)
    {
        U[q]  = new field4(p);
        U0[q] = new field4(p);
        R[q]  = new field4(p);
        Uc[q] = new field4(p);
    }
    B = new field4(p);
    phi = new field4(p);
    fallback_cells=0;
    plagoon=nullptr;
    print_idx=0;
    restart_time=0.0;

    gamma=p->A800;
    recon=p->A801;
    flux=p->A802;
    casetype=p->A803;
    dir=p->A804;
    x0=p->A805;
    curvature=p->A806;
    nuweno=p->A807;
    noslip=p->A808;
    bc[0]=p->A810_xm; bc[1]=p->A810_xp;
    bc[2]=p->A810_ym; bc[3]=p->A810_yp;
    bc[4]=p->A810_zm; bc[5]=p->A810_zp;
    rhoL=p->A811_rho; uL=p->A811_u; pL=p->A811_p;
    rhoR=p->A812_rho; uR=p->A812_u; pR=p->A812_p;
    zwall_from=p->A814;
    mu=p->A820; prandtl=p->A821; sutherland=p->A822;
    shearU=p->A823_U; shearH=p->A823_H;
    chanU=p->A824_U; chanH=p->A824_H; chanAng=p->A824_ang; chanN0=p->A824_n0;
    blastE=p->A825_E; blastR=p->A825_r0; blastX=p->A825_x; blastZ=p->A825_z;
    blastP0=p->A825_p0; blastRho0=p->A825_rho0;
    sponge_in=p->A826_in; sponge_out=p->A826_out; sponge_sigma=p->A826_sigma;
    sponge_mode=p->A827;
    top_x0=p->A828_x0;
    corner_x=p->A829_x; corner_z=p->A829_z; corner_n=p->A829_n;
    usenormal=p->A830;
    nrm[0]=p->A830_nx; nrm[1]=p->A830_ny; nrm[2]=p->A830_nz;
    cfl=p->N47;
    gridsolids=0;

    read_input();

    if(usenormal)
    {
        const double l=sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]);
        if(l>0.0)
        for(int q=0; q<3; ++q)
        nrm[q]/=l;
    }

    pack_state(p->A813_rho,p->A813_V,p->A813_ang,p->A813_p,Qfs);
    pack_state(p->A828_rho,p->A828_V,p->A828_ang,p->A828_p,Qtop);
}

comp_euler::~comp_euler()
{
    for(int q=0; q<5; ++q)
    {
        delete U[q];
        delete U0[q];
        delete R[q];
        delete Uc[q];
    }
    delete B;
    delete phi;
    delete plagoon;
}

void comp_euler::pack_state(double r, double V, double ang, double pr, double *Q) const
{
    const double a=ang*PI/180.0;
    Q[0]=r;
    Q[1]=r*V*cos(a);
    Q[2]=0.0;
    Q[3]=r*V*sin(a);
    Q[4]=pr/(gamma-1.0)+0.5*r*V*V;
}

void comp_euler::ini()
{
    p->count=0;
    p->printcount=0;
    p->printtime=0.0;
    cfl=p->N47;

    if(p->mpirank==0)
    {
        cout<<"REEF3D::COMP: compressible "<<(mu>0.0 ? "Navier-Stokes, viscosity "+std::to_string(mu)+", Prandtl "+std::to_string(prandtl)+"," : "Euler,")<<" ideal gas gamma "<<gamma
            <<", "<<(recon==1 ? "WENO5-Z characteristic" : "first order")<<", "<<(flux==1 ? "HLLC" : "HLL")
            <<", SSP-RK3, CFL "<<cfl<<endl;
        cout<<"grid "<<p->gknox<<" x "<<p->gknoy<<" x "<<p->gknoz<<endl;

        double QL[5]={1.0,0.3,0.1,-0.2,2.5}, QR[5]={0.5,-0.1,0.0,0.1,1.1}, Lm[5][5], Rm[5][5], err=0.0;
        eigen(QL,QR,Lm,Rm);
        for(int a=0; a<5; ++a)
        for(int b=0; b<5; ++b)
        {
            double s=0.0;
            for(int c=0; c<5; ++c)
            s+=Lm[a][c]*Rm[c][b];
            err=max(err,fabs(s-(a==b?1.0:0.0)));
        }
        cout<<"eigenvector check |L R - I| = "<<err<<endl;
    }

    solid_ini();
    initialize();
    boundary(U);

    if((!prints.empty() && prints[0]<=0.0) || (prints.empty() && p->P30>0.0))
    {
        print(p->printcount++);
        if(!prints.empty() && prints[0]<=0.0)
        ++print_idx;
        if(prints.empty() && p->P30>0.0)
        p->printtime+=p->P30;
    }
}

int comp_euler::advance()
{
    p->dt=timestep();
    double target=p->N41;
    if(print_idx<(int)prints.size())
    target=min(target,prints[print_idx]);
    p->dt=min(p->dt,max(target-p->simtime,1.0e-15));

    for(int q=0; q<5; ++q)
    ILOOP
    JLOOP
    KLOOP
    (*U0[q])(i,j,k)=(*U[q])(i,j,k);

    const double a[3]={0.0,0.75,1.0/3.0};
    for(int stage=0; stage<3; ++stage)
    {
        boundary(U);

        IMALOOP
        JMALOOP
        KMALOOP
        (*B)(i,j,k)=0.0;

        for(int iter=0; iter<6; ++iter)
        {
            rhs(U);

            for(int q=0; q<5; ++q)
            ILOOP
            JLOOP
            KLOOP
            (*Uc[q])(i,j,k) = solid(i,j,k) ? (*U[q])(i,j,k)
                               : a[stage]*(*U0[q])(i,j,k) + (1.0-a[stage])*((*U[q])(i,j,k) + p->dt*(*R[q])(i,j,k));

            if(recon==0 || check_stage(Uc)==0)
            break;

            pgc->gcparax(p,*B,4);
            ILOOP
            JLOOP
            KLOOP
            if((*B)(i,j,k)>0.0)
            ++fallback_cells;
        }

        for(int q=0; q<5; ++q)
        ILOOP
        JLOOP
        KLOOP
        (*U[q])(i,j,k) = (*Uc[q])(i,j,k);

        corner_fix(U);
    }

    double bad=0.0;
    ILOOP
    JLOOP
    KLOOP
    {
        double Q[5], r,u,v,w,pr;
        for(int q=0; q<5; ++q)
        Q[q]=(*U[q])(i,j,k);
        prim(Q,r,u,v,w,pr);
        if(!solid(i,j,k) && (!(r>0.0) || !(pr>0.0)))
        bad=1.0;
    }
    if(pgc->globalmax(bad)>0.0)
    {
        if(p->mpirank==0)
        cout<<"REEF3D::COMP: non-positive density or pressure at t = "<<p->simtime+p->dt<<endl;
        boundary(U);
        print(p->printcount++);
        return 1;
    }

    return 0;
}

void comp_euler::printout()
{
    int do_print=0;

    if(!prints.empty())
    {
        while(print_idx<(int)prints.size() && p->simtime>=prints[print_idx]-1.0e-12)
        {
            do_print=1;
            ++print_idx;
        }
    }
    else
    {
        if(p->P20>0 && p->P30<0.0 && p->count%p->P20==0)
        do_print=1;

        if(p->P30>0.0 && p->simtime>p->printtime)
        {
            do_print=1;
            p->printtime+=p->P30;
        }
    }

    if(do_print)
    print(p->printcount++);
}

void comp_euler::finish()
{
    if(plagoon)
    plagoon->finish(p,pgc);

    if(p->mpirank==0)
    cout<<"REEF3D::COMP finished: t "<<p->simtime<<" steps "<<p->count<<endl;
    double fb=pgc->globalsum(double(fallback_cells));
    if(p->mpirank==0)
    cout<<"REEF3D::COMP first-order fallback: "<<fb<<" cell-stages"<<endl;
}
