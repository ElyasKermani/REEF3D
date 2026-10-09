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

double comp_euler::muT(double T) const
{
    if(sutherland<=0.0)
    return mu;
    double r,u,v,w,pinf;
    prim(Qfs,r,u,v,w,pinf);
    const double th=T/(pinf/r);
    return mu*pow(th,1.5)*(1.0+sutherland)/(th+sutherland);
}

// Navier-Stokes fluxes at the faces, second-order central: stresses with Stokes' hypothesis and
// Fourier heat flux, k = mu cp / Pr with cp = gamma/(gamma-1) (R = 1); added to R as
// -(F_inv - F_visc). Cross derivatives average the central differences of the two face cells.
void comp_euler::viscous(field4 **V)
{
    if(mu<=0.0)
    return;

    const int gk[3]={p->gknox,p->gknoy,p->gknoz};
    const int nmax[3]={p->knox,p->knoy,p->knoz};

    // u v w T of every cell including the ghost layers, computed once
    const int jlo = p->gknoy>1 ? -3 : 0, jhi = p->gknoy>1 ? p->knoy+3 : p->knoy;
    const int NI=p->knox+6, NJ=jhi-jlo, NK=p->knoz+6;
    visc_prim.assign(size_t(NI)*NJ*NK*4,0.0);
    auto at=[&](int a, int b, int c) { return ((size_t(a+3)*NJ + size_t(b-jlo))*NK + size_t(c+3))*4; };
    for(int a=-3; a<p->knox+3; ++a)
    for(int b=jlo; b<jhi; ++b)
    for(int c=-3; c<p->knoz+3; ++c)
    {
        double Q[5], r,u,v,w,pr;
        for(int n=0; n<5; ++n) Q[n]=(*V[n])(a,b,c);
        prim(Q,r,u,v,w,pr);
        const size_t o=at(a,b,c);
        visc_prim[o]=u; visc_prim[o+1]=v; visc_prim[o+2]=w; visc_prim[o+3]=pr/r;
    }
    auto pr5=[&](int a, int b, int c, double *q)   // u v w T
    {
        const size_t o=at(a,b,c);
        q[0]=visc_prim[o]; q[1]=visc_prim[o+1]; q[2]=visc_prim[o+2]; q[3]=visc_prim[o+3];
    };
    auto hgt=[&](int d, int a, int b, int c)
    {
        return d==0 ? p->DXN[a+marge] : (d==1 ? p->DYN[b+marge] : p->DZN[c+marge]);
    };

    for(int d=0; d<3; ++d)
    {
        if(gk[d]<2)
        continue;

        const int a1 = d==0 ? p->knoy : p->knox;
        const int a2 = d==2 ? p->knoy : p->knoz;

        for(int q1=0; q1<a1; ++q1)
        for(int q2=0; q2<a2; ++q2)
        for(int f=0; f<=nmax[d]; ++f)
        {
            int L[3], Rr[3];
            if(d==0) { L[0]=f-1; Rr[0]=f; L[1]=Rr[1]=q1; L[2]=Rr[2]=q2; }
            else if(d==1) { L[1]=f-1; Rr[1]=f; L[0]=Rr[0]=q1; L[2]=Rr[2]=q2; }
            else { L[2]=f-1; Rr[2]=f; L[0]=Rr[0]=q1; L[1]=Rr[1]=q2; }

            // skip faces inside bodies
            if(solid(L[0],L[1],L[2]) && solid(Rr[0],Rr[1],Rr[2]))
            continue;

            double qL[4], qR[4];
            pr5(L[0],L[1],L[2],qL);
            pr5(Rr[0],Rr[1],Rr[2],qR);

            // gradients at the face: g[j][m], derivative of (u,v,w,T)[m] along x_j
            double g[3][4];
            const double hn=0.5*(hgt(d,L[0],L[1],L[2])+hgt(d,Rr[0],Rr[1],Rr[2]));
            for(int m=0; m<4; ++m)
            g[d][m]=(qR[m]-qL[m])/hn;

            for(int t=0; t<3; ++t)
            {
                if(t==d)
                continue;
                if(gk[t]<2)
                {
                    for(int m=0; m<4; ++m) g[t][m]=0.0;
                    continue;
                }
                int Lp[3]={L[0],L[1],L[2]}, Lm[3]={L[0],L[1],L[2]}, Rp[3]={Rr[0],Rr[1],Rr[2]}, Rm[3]={Rr[0],Rr[1],Rr[2]};
                Lp[t]++; Lm[t]--; Rp[t]++; Rm[t]--;
                double a[4],b[4],c[4],e[4];
                pr5(Lp[0],Lp[1],Lp[2],a); pr5(Lm[0],Lm[1],Lm[2],b); pr5(Rp[0],Rp[1],Rp[2],c); pr5(Rm[0],Rm[1],Rm[2],e);
                const double ht=hgt(t,L[0],L[1],L[2]);
                for(int m=0; m<4; ++m)
                g[t][m]=(a[m]-b[m]+c[m]-e[m])/(4.0*ht);
            }

            const double divu=g[0][0]+g[1][1]+g[2][2];
            const double muf=0.5*(muT(qL[3])+muT(qR[3]));
            const double kcond=muf*gamma/((gamma-1.0)*prandtl);
            double uf[3]={0.5*(qL[0]+qR[0]),0.5*(qL[1]+qR[1]),0.5*(qL[2]+qR[2])};
            double tau[3];
            for(int j=0; j<3; ++j)
            tau[j]=muf*(g[j][d]+g[d][j]) - (j==d ? 2.0/3.0*muf*divu : 0.0);

            double Fv[5];
            Fv[0]=0.0;
            Fv[1]=tau[0]; Fv[2]=tau[1]; Fv[3]=tau[2];
            Fv[4]=uf[0]*tau[0]+uf[1]*tau[1]+uf[2]*tau[2] + kcond*g[d][3];

            for(int side=0; side<2; ++side)
            {
                const int *C = side==0 ? L : Rr;
                if(C[d]<0 || C[d]>=nmax[d])
                continue;
                const double h=hgt(d,C[0],C[1],C[2]);
                const double sgn = side==0 ? 1.0 : -1.0;
                for(int n=0; n<5; ++n)
                (*R[n])(C[0],C[1],C[2]) += sgn*Fv[n]/h;
            }
        }
    }
}

// Woodward and Colella (1984): the flow round a sharp convex corner is nearly steady and isentropic, so the cells
// of the first row behind the corner take the entropy and total enthalpy of the cell just upstream of the corner
// (below and in front of it); their pressure and flow direction are kept. Removes the spurious entropy layer
// that the corner singularity otherwise produces and that grows with resolution.
