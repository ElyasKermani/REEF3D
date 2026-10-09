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
#include"geo_mesh.h"
#include<vector>
#include<cstring>
#include<cmath>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<sys/stat.h>

using namespace std;

void comp_euler::boundary(field4 **V)
{
    for(int n=0; n<5; ++n)
    pgc->gcparax(p,*V[n],4);

    const int m=3;
    const bool lo[3]={p->origin_i==0, p->origin_j==0, p->origin_k==0};
    const bool hi[3]={p->origin_i+p->knox==p->gknox, p->origin_j+p->knoy==p->gknoy, p->origin_k+p->knoz==p->gknoz};
    const int nmax[3]={p->knox,p->knoy,p->knoz};

    for(int d=0; d<3; ++d)
    {
        if(nmax[d]<1)
        continue;

        for(int side=0; side<2; ++side)
        {
            if((side==0 && !lo[d]) || (side==1 && !hi[d]))
            continue;

            const int btype=bc[2*d+side];
            // tangential ranges: directions treated before this one (and active) include their ghost layers,
            // so edge and corner ghosts are filled by the successive conditions
            const int ext_i = (d>0) ? m : 0;
            const int ext_j = (d>1 && p->gknoy>1) ? m : 0;
            const int q1lo = (d==0) ? 0 : -ext_i, q1hi = (d==0) ? p->knoy : p->knox+ext_i;
            const int q2lo = (d==2) ? -ext_j : 0,  q2hi = (d==2) ? p->knoy+ext_j : p->knoz;

            for(int q1=q1lo; q1<q1hi; ++q1)
            for(int q2=q2lo; q2<q2hi; ++q2)
            for(int g=1; g<=m; ++g)
            {
                int type=btype;
                if(d==2 && side==0 && p->XP[q1+marge]<zwall_from)
                {
                    if(type==1) type=2;        // slip wall downstream: free stream upstream (wedge frame)
                    else if(type==3) type=1;   // no-slip plate downstream: symmetry plane upstream
                }

                int ghost = side==0 ? -g : nmax[d]-1+g;
                int src = side==0 ? ((type==1||type==3) ? g-1 : 0) : ((type==1||type==3) ? nmax[d]-g : nmax[d]-1);

                int gi,gj,gk,si,sj,sk;
                if(d==0) { gi=ghost; si=src; gj=sj=q1; gk=sk=q2; }
                else if(d==1) { gj=ghost; sj=src; gi=si=q1; gk=sk=q2; }
                else { gk=ghost; sk=src; gi=si=q1; gj=sj=q2; }

                if(type==2)
                {
                    const bool shocked = d==2 && side==1 && p->XP[q1+marge]>=top_x0;
                    for(int n=0; n<5; ++n)
                    (*V[n])(gi,gj,gk) = shocked ? Qtop[n] : Qfs[n];
                    continue;
                }

                if(type==4 || type==5)
                {
                    // subsonic inflow (4): density and velocity of the free stream, pressure from inside;
                    // subsonic outflow (5): free-stream pressure, density and velocity from inside
                    double Qi[5], Qf[5], ri,ui,vi,wi,pi_,rf,uf,vf,wf,pf;
                    for(int n=0; n<5; ++n) { Qi[n]=(*V[n])(si,sj,sk); Qf[n]=Qfs[n]; }
                    prim(Qi,ri,ui,vi,wi,pi_);
                    prim(Qf,rf,uf,vf,wf,pf);
                    double r,u,v,w,pr;
                    if(type==4) { r=rf; u=uf; v=vf; w=wf; pr=pi_; }
                    else        { r=ri; u=ui; v=vi; w=wi; pr=pf; }
                    (*V[0])(gi,gj,gk)=r; (*V[1])(gi,gj,gk)=r*u; (*V[2])(gi,gj,gk)=r*v; (*V[3])(gi,gj,gk)=r*w;
                    (*V[4])(gi,gj,gk)=pr/(gamma-1.0)+0.5*r*(u*u+v*v+w*w);
                    continue;
                }

                for(int n=0; n<5; ++n)
                (*V[n])(gi,gj,gk) = (*V[n])(si,sj,sk);

                if(type==1 || type==3)
                (*V[1+d])(gi,gj,gk) = -(*V[1+d])(si,sj,sk);

                if(type==3)
                for(int m=1; m<=3; ++m)
                if(m!=1+d)
                (*V[m])(gi,gj,gk) = -(*V[m])(si,sj,sk);
            }
        }
    }

    if(!polys.empty() || gridsolids)
    {
        solid_fill(V);
        for(int n=0; n<5; ++n)
        pgc->gcparax(p,*V[n],4);
    }

}

// Roe-averaged eigenvectors in the local frame (normal velocity in slot 1)
void comp_euler::corners(field4 **V)
{
    for(int ii=-3; ii<p->knox+3; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=-3; kk<p->knoz+3; ++kk)
    {
        if((ii>=0 && ii<p->knox) || (kk>=0 && kk<p->knoz))
        continue;
        const int ck=max(0,min(p->knoz-1,kk));
        for(int n=0; n<5; ++n)
        (*V[n])(ii,jj,kk)=(*V[n])(ii,jj,ck);
    }
}

void comp_euler::corner_fix(field4 **V)
{
    if(corner_n<=0)
    return;

    auto find=[&](double x, double z, int &ci, int &ck)
    {
        ci=-1; ck=-1;
        for(int ii=0; ii<p->knox; ++ii) if(fabs(p->XP[ii+marge]-x)<=0.5*p->DXN[ii+marge]+1.0e-12) ci=ii;
        for(int kk=0; kk<p->knoz; ++kk) if(fabs(p->ZP[kk+marge]-z)<=0.5*p->DZN[kk+marge]+1.0e-12) ck=kk;
        return ci>=0 && ck>=0;
    };

    // reference: the fluid cell below and in front of the corner (single rank holds it)
    const double h=p->DXM;
    int ri,rk;
    double sref=0.0, Href=0.0, have=0.0;
    if(find(corner_x-0.5*h,corner_z-0.5*h,ri,rk))
    {
        double Q[5], r,u,v,w,pr;
        for(int n=0; n<5; ++n) Q[n]=(*V[n])(ri,0,rk);
        prim(Q,r,u,v,w,pr);
        sref=pr/pow(r,gamma);
        Href=gamma/(gamma-1.0)*pr/r+0.5*(u*u+v*v+w*w);
        have=1.0;
    }
    have=pgc->globalmax(have);
    sref=pgc->globalmax(sref);
    Href=pgc->globalmax(Href);
    if(have<=0.0)
    return;

    for(int q=0; q<corner_n; ++q)
    {
        int ci,ck;
        if(!find(corner_x+(q+0.5)*h,corner_z+0.5*h,ci,ck))
        continue;
        for(int jj=0; jj<p->knoy; ++jj)
        {
            double Q[5], r,u,v,w,pr;
            for(int n=0; n<5; ++n) Q[n]=(*V[n])(ci,jj,ck);
            prim(Q,r,u,v,w,pr);
            const double rn=pow(pr/sref,1.0/gamma);
            const double q2=max(0.0,2.0*(Href-gamma/(gamma-1.0)*pr/rn));
            const double sp=sqrt(u*u+v*v+w*w), sc = sp>1.0e-12 ? sqrt(q2)/sp : 0.0;
            u*=sc; v*=sc; w*=sc;
            (*V[0])(ci,jj,ck)=rn; (*V[1])(ci,jj,ck)=rn*u; (*V[2])(ci,jj,ck)=rn*v; (*V[3])(ci,jj,ck)=rn*w;
            (*V[4])(ci,jj,ck)=pr/(gamma-1.0)+0.5*rn*(u*u+v*v+w*w);
        }
    }
}

