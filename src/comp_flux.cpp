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

void comp_euler::prim(const double *Q, double &rho, double &u, double &v, double &w, double &pr) const
{
    rho=Q[0];
    u=Q[1]/rho;
    v=Q[2]/rho;
    w=Q[3]/rho;
    pr=(gamma-1.0)*(Q[4]-0.5*rho*(u*u+v*v+w*w));
}

void comp_euler::eigen(const double *QL, const double *QR, double (*Lm)[5], double (*Rm)[5]) const
{
    double rl,ul,vl,wl,pl,rr,ur,vr,wr,pr;
    prim(QL,rl,ul,vl,wl,pl);
    prim(QR,rr,ur,vr,wr,pr);

    const double sl=sqrt(rl), sr=sqrt(rr), f=1.0/(sl+sr);
    const double u=(sl*ul+sr*ur)*f, v=(sl*vl+sr*vr)*f, w=(sl*wl+sr*wr)*f;
    const double H=(sl*(QL[4]+pl)/rl + sr*(QR[4]+pr)/rr)*f;
    const double q2=u*u+v*v+w*w;
    const double c=sqrt(max((gamma-1.0)*(H-0.5*q2),1.0e-14));

    const double Rv[5][5]={{1.0,     1.0,    0.0, 0.0, 1.0},
                           {u-c,     u,      0.0, 0.0, u+c},
                           {v,       v,      1.0, 0.0, v},
                           {w,       w,      0.0, 1.0, w},
                           {H-u*c,   0.5*q2, v,   w,   H+u*c}};

    const double b1=(gamma-1.0)/(c*c), b2=0.5*b1*q2;
    const double Lv[5][5]={{0.5*(b2+u/c), -0.5*(b1*u+1.0/c), -0.5*b1*v, -0.5*b1*w, 0.5*b1},
                           {1.0-b2,        b1*u,               b1*v,      b1*w,     -b1},
                           {-v,            0.0,                1.0,       0.0,       0.0},
                           {-w,            0.0,                0.0,       1.0,       0.0},
                           {0.5*(b2-u/c), -0.5*(b1*u-1.0/c), -0.5*b1*v, -0.5*b1*w, 0.5*b1}};

    for(int a=0; a<5; ++a)
    for(int b=0; b<5; ++b)
    {
        Rm[a][b]=Rv[a][b];
        Lm[a][b]=Lv[a][b];
    }
}

// WENO5-Z value at the face between c and d, upwind cell c (Borges et al. 2008)
void comp_euler::hllc(const double *QL, const double *QR, double *F)
{
    double rl,ul,vl,wl,pl,rr,ur,vr,wr,pr;
    prim(QL,rl,ul,vl,wl,pl);
    prim(QR,rr,ur,vr,wr,pr);
    const double cl=sqrt(gamma*pl/rl), cr=sqrt(gamma*pr/rr);

    const double sl=sqrt(rl), sr=sqrt(rr), f=1.0/(sl+sr);
    const double ut=(sl*ul+sr*ur)*f, vt=(sl*vl+sr*vr)*f, wt=(sl*wl+sr*wr)*f;
    const double Ht=(sl*(QL[4]+pl)/rl + sr*(QR[4]+pr)/rr)*f;
    const double ct=sqrt(max((gamma-1.0)*(Ht-0.5*(ut*ut+vt*vt+wt*wt)),1.0e-14));

    const double SL=min(ul-cl,ut-ct), SR=max(ur+cr,ut+ct);
    const double SM=(pr-pl + rl*ul*(SL-ul) - rr*ur*(SR-ur))/(rl*(SL-ul) - rr*(SR-ur));

    const double FL[5]={rl*ul, rl*ul*ul+pl, rl*ul*vl, rl*ul*wl, ul*(QL[4]+pl)};
    const double FR[5]={rr*ur, rr*ur*ur+pr, rr*ur*vr, rr*ur*wr, ur*(QR[4]+pr)};

    if(SL>=0.0)
    {
        for(int n=0; n<5; ++n) F[n]=FL[n];
    }
    else if(SR<=0.0)
    {
        for(int n=0; n<5; ++n) F[n]=FR[n];
    }
    else if(SM>=0.0)
    {
        const double fac=rl*(SL-ul)/(SL-SM);
        const double Qs[5]={fac, fac*SM, fac*vl, fac*wl, fac*(QL[4]/rl + (SM-ul)*(SM + pl/(rl*(SL-ul))))};
        for(int n=0; n<5; ++n) F[n]=FL[n]+SL*(Qs[n]-QL[n]);
    }
    else
    {
        const double fac=rr*(SR-ur)/(SR-SM);
        const double Qs[5]={fac, fac*SM, fac*vr, fac*wr, fac*(QR[4]/rr + (SM-ur)*(SM + pr/(rr*(SR-ur))))};
        for(int n=0; n<5; ++n) F[n]=FR[n]+SR*(Qs[n]-QR[n]);
    }
}

void comp_euler::hll(const double *QL, const double *QR, double *F)
{
    double rl,ul,vl,wl,pl,rr,ur,vr,wr,pr;
    prim(QL,rl,ul,vl,wl,pl);
    prim(QR,rr,ur,vr,wr,pr);
    const double cl=sqrt(gamma*pl/rl), cr=sqrt(gamma*pr/rr);

    const double sl=sqrt(rl), sr=sqrt(rr), f=1.0/(sl+sr);
    const double ut=(sl*ul+sr*ur)*f, vt=(sl*vl+sr*vr)*f, wt=(sl*wl+sr*wr)*f;
    const double Ht=(sl*(QL[4]+pl)/rl + sr*(QR[4]+pr)/rr)*f;
    const double ct=sqrt(max((gamma-1.0)*(Ht-0.5*(ut*ut+vt*vt+wt*wt)),1.0e-14));
    const double SL=min(ul-cl,ut-ct), SR=max(ur+cr,ut+ct);

    const double FL[5]={rl*ul, rl*ul*ul+pl, rl*ul*vl, rl*ul*wl, ul*(QL[4]+pl)};
    const double FR[5]={rr*ur, rr*ur*ur+pr, rr*ur*vr, rr*ur*wr, ur*(QR[4]+pr)};

    for(int n=0; n<5; ++n)
    {
        if(SL>=0.0) F[n]=FL[n];
        else if(SR<=0.0) F[n]=FR[n];
        else F[n]=(SR*FL[n]-SL*FR[n]+SL*SR*(QR[n]-QL[n]))/(SR-SL);
    }
}

// Qs: six cells along direction d, the face lies between Qs[2] and Qs[3] (global ordering)
void comp_euler::face_flux(const double *fo, const double (*Qs)[5], int d, double *F)
{
    face_flux_h(fo,Qs,d,F,nullptr);
}

void comp_euler::face_flux_h(const double *fo, const double (*Qs)[5], int d, double *F, const double *hw)
{
    const int idx[5]={0, 1+d, 1+(d+1)%3, 1+(d+2)%3, 4};
    double W[6][5], QL[5], QR[5], Fl[5];

    for(int s=0; s<6; ++s)
    for(int n=0; n<5; ++n)
    W[s][n]=Qs[s][idx[n]];

    for(int n=0; n<5; ++n)
    {
        QL[n]=W[2][n];
        QR[n]=W[3][n];
    }

    if(recon==1 && !(fo && *fo>0.0))
    {
        double Lm[5][5], Rm[5][5], C[6][5], CL[5], CR[5], qL[5], qR[5];
        eigen(W[2],W[3],Lm,Rm);

        for(int s=0; s<6; ++s)
        for(int a=0; a<5; ++a)
        {
            C[s][a]=0.0;
            for(int b=0; b<5; ++b)
            C[s][a]+=Lm[a][b]*W[s][b];
        }

        bool uniform=true;
        if(hw && nuweno==1)
        for(int q=0; q<6; ++q)
        if(fabs(hw[q]-hw[2])>1.0e-10*hw[2]) uniform=false;

        nucoef kl, kr;
        if(!uniform)
        {
            const double hl[5]={hw[0],hw[1],hw[2],hw[3],hw[4]}, hr[5]={hw[5],hw[4],hw[3],hw[2],hw[1]};
            nu_coefficients(hl,kl);
            nu_coefficients(hr,kr);
        }
        for(int a=0; a<5; ++a)
        {
            if(uniform)
            {
                CL[a]=weno5z(C[0][a],C[1][a],C[2][a],C[3][a],C[4][a]);
                CR[a]=weno5z(C[5][a],C[4][a],C[3][a],C[2][a],C[1][a]);
            }
            else
            {
                const double vl[5]={C[0][a],C[1][a],C[2][a],C[3][a],C[4][a]};
                const double vr[5]={C[5][a],C[4][a],C[3][a],C[2][a],C[1][a]};
                CL[a]=nu_weno(vl,kl);
                CR[a]=nu_weno(vr,kr);
            }
        }

        for(int a=0; a<5; ++a)
        {
            qL[a]=qR[a]=0.0;
            for(int b=0; b<5; ++b)
            {
                qL[a]+=Rm[a][b]*CL[b];
                qR[a]+=Rm[a][b]*CR[b];
            }
        }

        double r,u,v,w,pr;
        prim(qL,r,u,v,w,pr);
        if(r>0.0 && pr>0.0 && std::isfinite(pr))
        for(int n=0; n<5; ++n) QL[n]=qL[n];

        prim(qR,r,u,v,w,pr);
        if(r>0.0 && pr>0.0 && std::isfinite(pr))
        for(int n=0; n<5; ++n) QR[n]=qR[n];
    }

    if(flux==1)
    hllc(QL,QR,Fl);
    else
    hll(QL,QR,Fl);

    for(int n=0; n<5; ++n)
    F[idx[n]]=Fl[n];
}

void comp_euler::rhs(field4 **V)
{
    for(int n=0; n<5; ++n)
    ILOOP
    JLOOP
    KLOOP
    (*R[n])(i,j,k)=0.0;

    const int gk[3]={p->gknox,p->gknoy,p->gknoz};
    const int nmax[3]={p->knox,p->knoy,p->knoz};
    double Qs[6][5], F[5];

    for(int d=0; d<3; ++d)
    {
        if(gk[d]<2)
        continue;

        const int a1 = d==0 ? p->knoy : p->knox;
        const int a2 = d==2 ? p->knoy : p->knoz;

        for(int q1=0; q1<a1; ++q1)
        for(int q2=0; q2<a2; ++q2)
        for(int f=0; f<=nmax[d]; ++f)     // face between cells f-1 and f
        {
            for(int s=0; s<6; ++s)
            {
                int c=f-3+s, ci,cj,ck;
                if(d==0) { ci=c; cj=q1; ck=q2; }
                else if(d==1) { cj=c; ci=q1; ck=q2; }
                else { ck=c; ci=q1; cj=q2; }

                for(int n=0; n<5; ++n)
                Qs[s][n]=(*V[n])(ci,cj,ck);
            }

            double fo=0.0;
            {
                int cl=f-1, cr=f, li,lj,lk,ri,rj,rk;
                if(d==0) { li=cl; ri=cr; lj=rj=q1; lk=rk=q2; }
                else if(d==1) { lj=cl; rj=cr; li=ri=q1; lk=rk=q2; }
                else { lk=cl; rk=cr; li=ri=q1; lj=rj=q2; }
                fo=max((*B)(li,lj,lk),(*B)(ri,rj,rk));
            }
            double hw[6];
            for(int q=0; q<6; ++q)
            {
                const int c=f-3+q;
                hw[q] = d==0 ? p->DXN[c+marge] : (d==1 ? p->DYN[c+marge] : p->DZN[c+marge]);
            }
            face_flux_h(&fo,Qs,d,F,hw);

            for(int side=0; side<2; ++side)
            {
                int c = f-1+side;
                if(c<0 || c>=nmax[d])
                continue;

                int ci,cj,ck;
                if(d==0) { ci=c; cj=q1; ck=q2; }
                else if(d==1) { cj=c; ci=q1; ck=q2; }
                else { ck=c; ci=q1; cj=q2; }

                double h = d==0 ? p->DXN[ci+marge] : (d==1 ? p->DYN[cj+marge] : p->DZN[ck+marge]);
                double sgn = side==0 ? -1.0 : 1.0;

                for(int n=0; n<5; ++n)
                (*R[n])(ci,cj,ck) += sgn*F[n]/h;
            }
        }
    }
    viscous(V);

    // absorbing zones: relaxation towards the free stream, sigma (d/L)^2
    if(sponge_sigma>0.0)
    {
        double xmin=p->XN[marge], xmax=p->XN[p->knox+marge];
        xmin=pgc->globalmin(xmin); xmax=pgc->globalmax(xmax);
        ILOOP
        JLOOP
        KLOOP
        {
            const double x=p->XP[IP];
            double sg=0.0;
            if(x<sponge_in) sg=sponge_sigma*pow((sponge_in-x)/(sponge_in-xmin),2.0);
            if(x>sponge_out) sg=sponge_sigma*pow((x-sponge_out)/(xmax-sponge_out),2.0);
            if(sg<=0.0)
            continue;

            if(sponge_mode==1)
            {
                // energy only: towards the free-stream pressure with the local density and velocity
                double Q[5], r,u,v,w,pr, rf,uf,vf,wf,pf;
                for(int n=0; n<5; ++n) Q[n]=(*V[n])(i,j,k);
                prim(Q,r,u,v,w,pr);
                prim(Qfs,rf,uf,vf,wf,pf);
                (*R[4])(i,j,k) -= sg*(pr-pf)/(gamma-1.0);
            }
            else
            for(int n=0; n<5; ++n)
            (*R[n])(i,j,k) -= sg*((*V[n])(i,j,k)-Qfs[n]);
        }
    }
}

// ghost cells outside the domain in two directions (x-z): copied from the ghost along x,
// needed by the cross derivatives of the viscous stresses at boundary faces
