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

bool comp_euler::solid(int a, int b, int c) const
{
    return (!polys.empty() || gridsolids) && (*phi)(a,b,c)<0.0;
}

// a surface edge lying on the domain boundary is not a wall the fluid sees (a body cut by the boundary)
bool comp_euler::on_boundary(double x1, double z1, double x2, double z2) const
{
    const double e=1.0e-6*(dom[1]-dom[0]);
    return (x1<=dom[0]+e && x2<=dom[0]+e) || (x1>=dom[1]-e && x2>=dom[1]-e)
        || (z1<=dom[2]+e && z2<=dom[2]+e) || (z1>=dom[3]-e && z2>=dom[3]-e);
}

// signed distance to the union of the polygons (negative inside), and the closest surface point
double comp_euler::solid_distance(double x, double z, double &bx, double &bz) const
{
    double best=1.0e30, sgn=1.0;
    bx=x; bz=z;
    for(const std::vector<double> &v : polys)
    {
        const int nv=int(v.size()/2);
        bool inside=false;
        double dmin=1.0e30, cx=x, cz=z;
        for(int q=0; q<nv; ++q)
        {
            const double x1=v[2*q], z1=v[2*q+1], x2=v[2*((q+1)%nv)], z2=v[2*((q+1)%nv)+1];
            if(((z1>z)!=(z2>z)) && (x < (x2-x1)*(z-z1)/(z2-z1)+x1))
            inside=!inside;
            if(on_boundary(x1,z1,x2,z2))
            continue;
            const double ex=x2-x1, ez=z2-z1, L2=ex*ex+ez*ez;
            double t = L2>0.0 ? ((x-x1)*ex+(z-z1)*ez)/L2 : 0.0;
            t=max(0.0,min(1.0,t));
            const double px=x1+t*ex, pz=z1+t*ez, d=sqrt((x-px)*(x-px)+(z-pz)*(z-pz));
            if(d<dmin) { dmin=d; cx=px; cz=pz; }
        }
        // union: inside any body wins; outside all, the nearest surface
        if(inside)
        {
            if(sgn>0.0 || dmin<best) { best=dmin; bx=cx; bz=cz; }
            sgn=-1.0;
        }
        else if(sgn>0.0 && dmin<best) { best=dmin; bx=cx; bz=cz; }
    }
    return sgn*best;
}

double comp_euler::segment_distance(double x, double z, double &bx, double &bz) const
{
    double best=1.0e30;
    bx=x; bz=z;
    for(size_t q=0; q+3<segs.size(); q+=4)
    {
        const double x1=segs[q], z1=segs[q+1], ex=segs[q+2]-x1, ez=segs[q+3]-z1, L2=ex*ex+ez*ez;
        double t = L2>0.0 ? ((x-x1)*ex+(z-z1)*ez)/L2 : 0.0;
        t=max(0.0,min(1.0,t));
        const double px=x1+t*ex, pz=z1+t*ez, d=sqrt((x-px)*(x-px)+(z-pz)*(z-pz));
        if(d<best) { best=d; bx=px; bz=pz; }
    }
    return best;
}

double comp_euler::triangle_distance(double x, double y, double z, double &bx, double &by, double &bz) const
{
    // closest point on each triangle (Ericson 2005, Real-Time Collision Detection 5.1.5)
    double best=1.0e30;
    bx=x; by=y; bz=z;
    for(size_t t=0; t+8<tris.size(); t+=9)
    {
        const double ax=tris[t],ay=tris[t+1],az=tris[t+2], ux=tris[t+3],uy=tris[t+4],uz=tris[t+5], vx=tris[t+6],vy=tris[t+7],vz=tris[t+8];
        const double abx=ux-ax,aby=uy-ay,abz=uz-az, acx=vx-ax,acy=vy-ay,acz=vz-az, apx=x-ax,apy=y-ay,apz=z-az;
        double cx,cy,cz;
        const double d1=abx*apx+aby*apy+abz*apz, d2=acx*apx+acy*apy+acz*apz;
        const double bpx=x-ux,bpy=y-uy,bpz=z-uz, cpx=x-vx,cpy=y-vy,cpz=z-vz;
        const double d3=abx*bpx+aby*bpy+abz*bpz, d4=acx*bpx+acy*bpy+acz*bpz;
        const double d5=abx*cpx+aby*cpy+abz*cpz, d6=acx*cpx+acy*cpy+acz*cpz;
        const double vc=d1*d4-d3*d2, vb=d5*d2-d1*d6, va=d3*d6-d5*d4;
        if(d1<=0.0 && d2<=0.0) { cx=ax; cy=ay; cz=az; }
        else if(d3>=0.0 && d4<=d3) { cx=ux; cy=uy; cz=uz; }
        else if(vc<=0.0 && d1>=0.0 && d3<=0.0) { const double w=d1/(d1-d3); cx=ax+w*abx; cy=ay+w*aby; cz=az+w*abz; }
        else if(d6>=0.0 && d5<=d6) { cx=vx; cy=vy; cz=vz; }
        else if(vb<=0.0 && d2>=0.0 && d6<=0.0) { const double w=d2/(d2-d6); cx=ax+w*acx; cy=ay+w*acy; cz=az+w*acz; }
        else if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) { const double w=(d4-d3)/((d4-d3)+(d5-d6)); cx=ux+w*(vx-ux); cy=uy+w*(vy-uy); cz=uz+w*(vz-uz); }
        else { const double den=1.0/(va+vb+vc), w1=vb*den, w2=vc*den; cx=ax+abx*w1+acx*w2; cy=ay+aby*w1+acy*w2; cz=az+abz*w1+acz*w2; }
        const double dd=(x-cx)*(x-cx)+(y-cy)*(y-cy)+(z-cz)*(z-cz);
        if(dd<best) { best=dd; bx=cx; by=cy; bz=cz; }
    }
    return sqrt(best);
}

void comp_euler::solid_ini()
{
    dom[0]=pgc->globalmin(p->XN[marge]); dom[1]=pgc->globalmax(p->XN[p->knox+marge]);
    dom[2]=pgc->globalmin(p->ZN[marge]); dom[3]=pgc->globalmax(p->ZN[p->knoz+marge]);

    // solids of the grid file: the signed distance of the geometry core (lexer::grid_solids),
    // used when comp.txt defines no polygons
    gridsolids=0;
    if(polys.empty() && p->flag_solid!=nullptr)
    {
        double neg=0.0;
        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        for(int kk=0; kk<p->knoz; ++kk)
        if(p->flag_solid[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0.0)
        neg=1.0;
        gridsolids = pgc->globalmax(neg)>0.0;
    }

    for(int ii=-3; ii<p->knox+3; ++ii)
    for(int jj=(p->gknoy>1?-3:0); jj<(p->gknoy>1?p->knoy+3:p->knoy); ++jj)
    for(int kk=-3; kk<p->knoz+3; ++kk)
    {
        double bx,bz;
        if(gridsolids && !segs.empty())
        segment_distance(p->XP[ii+marge],p->ZP[kk+marge],bx,bz);
        else if(gridsolids)
        {
            const int ci=max(0,min(p->knox-1,ii)), cj=max(0,min(p->knoy-1,jj)), ck=max(0,min(p->knoz-1,kk));
            (*phi)(ii,jj,kk) = p->flag_solid[(ci-p->imin)*p->jmax*p->kmax + (cj-p->jmin)*p->kmax + ck-p->kmin];
        }
        else
        {
            // beyond a physical domain boundary the ghost cells take the level set of the nearest
            // interior cell, so a body ending at the boundary does not look like fluid behind it
            int ai=ii, ak=kk;
            if(p->origin_i==0 && ai<0) ai=0;
            if(p->origin_i+p->knox==p->gknox && ai>=p->knox) ai=p->knox-1;
            if(p->origin_k==0 && ak<0) ak=0;
            if(p->origin_k+p->knoz==p->gknoz && ak>=p->knoz) ak=p->knoz-1;
            (*phi)(ii,jj,kk) = polys.empty() ? 1.0 : solid_distance(p->XP[ai+marge],p->ZP[ak+marge],bx,bz);
        }
    }

    if(gridsolids)
    pgc->gcparax(p,*phi,4);   // the ghost layers at subdomain interfaces from the neighbours

    // 3D: the exact distance to the solid triangles within 5 cells of the surface, sign from the geometry core
    if(gridsolids && p->gknoy>1 && p->gridgeo!=nullptr)
    {
        const geo_mesh &G=*p->gridgeo;
        for(const geo_object &o : G.obj)
        {
            if(o.role!=geo_mesh::role_solid) continue;
            for(int t=o.ts; t<o.te; ++t)
            for(int v=0; v<3; ++v)
            { tris.push_back(G.tri_x[t][v]); tris.push_back(G.tri_y[t][v]); tris.push_back(G.tri_z[t][v]); }
        }
        const double band=5.0*p->DXM;
        for(int ii=-3; ii<p->knox+3; ++ii)
        for(int jj=-3; jj<p->knoy+3; ++jj)
        for(int kk=-3; kk<p->knoz+3; ++kk)
        {
            double &ph=(*phi)(ii,jj,kk), bx,by,bz;
            if(fabs(ph)<band && !tris.empty())
            ph = (ph<0.0 ? -1.0 : 1.0)*triangle_distance(p->XP[ii+marge],p->YP[jj+marge],p->ZP[kk+marge],bx,by,bz);
        }
        const double ntri=pgc->globalmax(double(tris.size()/9));
        if(p->mpirank==0)
        cout<<"REEF3D::COMP solids: "<<ntri<<" surface triangles, exact distance within 5 cells (3D)"<<endl;
    }

    // 2D: the exact distance near the surface from the solid triangles cut by the plane of the
    // cells (the sign stays that of the geometry core, whose distance is not exact near the surface)
    if(gridsolids && p->gknoy==1 && p->gridgeo!=nullptr)
    {
        const geo_mesh &G=*p->gridgeo;
        const double yc=p->YP[marge]+1.0e-9;
        for(const geo_object &o : G.obj)
        {
            if(o.role!=geo_mesh::role_solid)
            continue;
            for(int t=o.ts; t<o.te; ++t)
            {
                double pts[4]; int np=0;
                for(int e=0; e<3; ++e)
                {
                    const int a=e, b=(e+1)%3;
                    const double ya=G.tri_y[t][a]-yc, yb=G.tri_y[t][b]-yc;
                    if((ya<0.0)!=(yb<0.0) && np<4)
                    {
                        const double f=ya/(ya-yb);
                        pts[np++]=G.tri_x[t][a]+f*(G.tri_x[t][b]-G.tri_x[t][a]);
                        pts[np++]=G.tri_z[t][a]+f*(G.tri_z[t][b]-G.tri_z[t][a]);
                    }
                }
                if(np==4 && !on_boundary(pts[0],pts[1],pts[2],pts[3]))
                segs.insert(segs.end(),pts,pts+4);
            }
        }

        const double band=5.0*p->DXM;
        for(int ii=-3; ii<p->knox+3; ++ii)
        for(int kk=-3; kk<p->knoz+3; ++kk)
        {
            double &ph=(*phi)(ii,0,kk), bx,bz;
            if(fabs(ph)<band && !segs.empty())
            ph = (ph<0.0 ? -1.0 : 1.0)*segment_distance(p->XP[ii+marge],p->ZP[kk+marge],bx,bz);
        }

        const double nseg=pgc->globalmax(double(segs.size()/4));
        if(p->mpirank==0)
        cout<<"REEF3D::COMP solids: "<<nseg<<" surface segments in the plane of the cells, exact distance within 5 cells"<<endl;
    }

    if(polys.empty() && !gridsolids)
    return;

    long nib=0, nsolid=0;
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    for(int kk=0; kk<p->knoz; ++kk)
    {
        if((*phi)(ii,jj,kk)>=0.0)
        continue;
        ++nsolid;

        // within reach of a fluid stencil (3 cells along x or z)?
        bool near=false;
        for(int s=-3; s<=3 && !near; ++s)
        {
            if(ii+s>=-3 && ii+s<p->knox+3 && (*phi)(ii+s,jj,kk)>=0.0) near=true;
            if(kk+s>=-3 && kk+s<p->knoz+3 && (*phi)(ii,jj,kk+s)>=0.0) near=true;
            if(p->gknoy>1 && jj+s>=-3 && jj+s<p->knoy+3 && (*phi)(ii,jj+s,kk)>=0.0) near=true;
        }

        if(!near)
        {
            deep.push_back(ii); deep.push_back(jj); deep.push_back(kk);
            continue;
        }

        double bx,bz, by=p->YP[jj+marge];
        if(gridsolids && !tris.empty())
        triangle_distance(p->XP[ii+marge],p->YP[jj+marge],p->ZP[kk+marge],bx,by,bz);
        else if(gridsolids && !segs.empty())
        segment_distance(p->XP[ii+marge],p->ZP[kk+marge],bx,bz);
        else if(gridsolids)
        {
            // closest surface point from the signed distance and its gradient
            const double ph=(*phi)(ii,jj,kk);
            double gx=((*phi)(ii+1,jj,kk)-(*phi)(ii-1,jj,kk))/(p->XP[ii+1+marge]-p->XP[ii-1+marge]);
            double gz=((*phi)(ii,jj,kk+1)-(*phi)(ii,jj,kk-1))/(p->ZP[kk+1+marge]-p->ZP[kk-1+marge]);
            const double gn=sqrt(gx*gx+gz*gz);
            if(gn>0.0) { gx/=gn; gz/=gn; } else { gx=0.0; gz=1.0; }
            bx=p->XP[ii+marge]-ph*gx;
            bz=p->ZP[kk+marge]-ph*gz;
        }
        else
        solid_distance(p->XP[ii+marge],p->ZP[kk+marge],bx,bz);
        ibcell c;
        c.i=ii; c.j=jj; c.k=kk;
        // image point on the normal through the closest wall point, at least 1.5 cells from the wall so that
        // its bilinear stencil lies in the fluid; the ghost value is then the wall condition linearly
        // extrapolated through the image value (second order, Mittal et al. 2008)
        {
            const double gx=p->XP[ii+marge], gy=p->YP[jj+marge], gz=p->ZP[kk+marge];
            double nx=bx-gx, ny=by-gy, nz=bz-gz;
            const double dg=sqrt(nx*nx+ny*ny+nz*nz);
            if(dg>1.0e-12*p->DXM) { nx/=dg; ny/=dg; nz/=dg; }
            else
            {
                nx=((*phi)(ii+1,jj,kk)-(*phi)(ii-1,jj,kk)); nz=((*phi)(ii,jj,kk+1)-(*phi)(ii,jj,kk-1));
                ny = p->gknoy>1 ? ((*phi)(ii,jj+1,kk)-(*phi)(ii,jj-1,kk)) : 0.0;
                const double l=sqrt(nx*nx+ny*ny+nz*nz); if(l>0.0) { nx/=l; ny/=l; nz/=l; } else { nx=0.0; ny=0.0; nz=1.0; }
            }
            c.ny=ny;
            const double hmin=min(p->DXN[ii+marge],p->DZN[kk+marge]);
            // slip walls (inviscid) keep the mirror point: near a sharp leading edge the flow changes across
            // a few cells and a distant image point misplaced the wedge shock by 0.06; no-slip walls use the
            // distant image point, which halved to sixfold reduced the error of the inclined channel test
            c.dg=dg;
            c.di = noslip==1 ? max(dg,1.5*hmin) : dg;
            c.nx=nx; c.nz=nz;
            c.xi=bx+c.di*nx;
            c.yi=by+c.di*ny;
            c.zi=bz+c.di*nz;
            c.dn=c.dg+c.di;
        }

        // wall curvature from the level set at the ghost cell: kappa_G = div(grad phi/|grad phi|),
        // shifted to the wall, kappa_w = kappa_G/(1 + d kappa_G), d = |phi|; limited to 1/(2 dx)
        {
            const double hx=p->DXN[ii+marge], hz=p->DZN[kk+marge];
            const double fx=((*phi)(ii+1,jj,kk)-(*phi)(ii-1,jj,kk))/(2.0*hx), fz=((*phi)(ii,jj,kk+1)-(*phi)(ii,jj,kk-1))/(2.0*hz);
            const double fxx=((*phi)(ii+1,jj,kk)-2.0*(*phi)(ii,jj,kk)+(*phi)(ii-1,jj,kk))/(hx*hx);
            const double fzz=((*phi)(ii,jj,kk+1)-2.0*(*phi)(ii,jj,kk)+(*phi)(ii,jj,kk-1))/(hz*hz);
            const double fxz=((*phi)(ii+1,jj,kk+1)-(*phi)(ii+1,jj,kk-1)-(*phi)(ii-1,jj,kk+1)+(*phi)(ii-1,jj,kk-1))/(4.0*hx*hz);
            const double g2=fx*fx+fz*fz;
            double kg = g2>1.0e-12 ? (fxx*fz*fz - 2.0*fx*fz*fxz + fzz*fx*fx)/pow(g2,1.5) : 0.0;
            const double d=fabs((*phi)(ii,jj,kk));
            double kw = fabs(1.0+d*kg)>1.0e-6 ? kg/(1.0+d*kg) : 0.0;
            // the correction assumes a smooth wall resolved by the grid: at sharp corners (curvature at the
            // limit, or a change of more than 25 % across the ghost-image distance) the plain mirror is used
            const double kmax=0.5/min(hx,hz);
            c.kappa = (fabs(kw)<kmax && fabs(kw)*c.dn<=0.25) ? kw : 0.0;

            // 3D: the curvature correction is not implemented (normal curvature along the flow), plain mirror
            if(p->gknoy>1) c.kappa=0.0;

            // the stencil above reaches past a physical domain boundary: no reliable curvature there
            if((ii==0 && p->origin_i==0) || (ii==p->knox-1 && p->origin_i+p->knox==p->gknox)
            || (kk==0 && p->origin_k==0) || (kk==p->knoz-1 && p->origin_k+p->knoz==p->gknoz))
            c.kappa=0.0;
        }
        ibcells.push_back(c);
        ++nib;
    }

    const double tib=pgc->globalsum(double(nib)), tsolid=pgc->globalsum(double(nsolid));
    if(p->mpirank==0)
    cout<<"REEF3D::COMP wall ghost cells: "<<(curvature==1 ? "curvature-corrected symmetry" : "plain mirror")<<endl;
    if(p->mpirank==0)
    cout<<"REEF3D::COMP solids: "<<(gridsolids ? "grid file (geometry core)" : std::to_string(polys.size())+" polygon bodies")<<", "<<tsolid<<" solid cells, "<<tib<<" immersed-boundary ghost cells"<<endl;
}

// ghost cells inside the bodies: the state at the image point (mirrored across the surface),
// with the normal velocity reversed (slip wall); bilinear interpolation over the fluid cells only
void comp_euler::solid_fill(field4 **V)
{
    if(polys.empty() && !gridsolids)
    return;

    const double x0=p->XP[marge], y0=p->YP[marge], z0=p->ZP[marge];
    const double dx=p->DXN[marge], dy=p->DYN[marge], dz=p->DZN[marge];
    const bool td = p->gknoy>1;

    for(const ibcell &c : ibcells)
    {
        const double fx=(c.xi-x0)/dx, fz=(c.zi-z0)/dz, fy = td ? (c.yi-y0)/dy : double(c.j);
        int i0=int(floor(fx)), k0=int(floor(fz)), j0 = td ? int(floor(fy)) : c.j;
        const double ax=fx-i0, az=fz-k0, ay = td ? fy-j0 : 0.0;
        double W=0.0, prim5[5]={0,0,0,0,0};

        for(int a=0; a<2; ++a)
        for(int e=0; e<(td?2:1); ++e)
        for(int b=0; b<2; ++b)
        {
            int ii=max(-3,min(p->knox+2,i0+a)), kk=max(-3,min(p->knoz+2,k0+b));
            int jj = td ? max(-3,min(p->knoy+2,j0+e)) : c.j;
            if((*phi)(ii,jj,kk)<0.0)
            continue;
            const double w=(a?ax:1.0-ax)*(b?az:1.0-az)*(td ? (e?ay:1.0-ay) : 1.0)+1.0e-12;
            double Q[5], r,u,v,ww,pr;
            for(int n=0; n<5; ++n) Q[n]=(*V[n])(ii,jj,kk);
            prim(Q,r,u,v,ww,pr);
            prim5[0]+=w*r; prim5[1]+=w*u; prim5[2]+=w*v; prim5[3]+=w*ww; prim5[4]+=w*pr;
            W+=w;
        }

        if(W<=0.0)
        {
            // no fluid cell around the image point: the nearest fluid neighbour within two cells
            double dbest=1.0e30;
            for(int a=-2; a<=3; ++a)
            for(int e=(td?-2:0); e<=(td?3:0); ++e)
            for(int b=-2; b<=3; ++b)
            {
                int ii=max(-3,min(p->knox+2,i0+a)), kk=max(-3,min(p->knoz+2,k0+b));
                int jj = td ? max(-3,min(p->knoy+2,j0+e)) : c.j;
                if((*phi)(ii,jj,kk)<0.0) continue;
                const double d=fabs(a-ax)+fabs(b-az)+(td?fabs(e-ay):0.0);
                if(d<dbest)
                {
                    dbest=d;
                    double Q[5]; for(int n=0; n<5; ++n) Q[n]=(*V[n])(ii,jj,kk);
                    prim(Q,prim5[0],prim5[1],prim5[2],prim5[3],prim5[4]);
                    W=1.0;
                }
            }
            if(W<=0.0) continue;
        }
        else
        for(int n=0; n<5; ++n) prim5[n]/=W;

        double r=prim5[0], u=prim5[1], v=prim5[2], w=prim5[3], pr=prim5[4];
        const double un=u*c.nx+v*c.ny+w*c.nz;
        const double f = c.di>0.0 ? c.dg/c.di : 1.0;   // linear extrapolation through the wall value

        if(noslip==1)
        {
            u=-f*u; v=-f*v; w=-f*w;
        }
        else if(curvature==1 && c.kappa!=0.0)
        {
            // curvature-corrected symmetry (Dadone & Grossman): normal momentum dp/dn = rho ut^2 kappa
            // across the ghost-image distance, constant entropy and total enthalpy along the normal
            const double tx=u-un*c.nx, tz=w-un*c.nz, ut2=tx*tx+v*v+tz*tz;
            double pg = pr - r*ut2*c.kappa*c.dn;
            pg = max(0.5*pr,min(1.5*pr,pg));
            const double rg = r*pow(pg/pr,1.0/gamma);
            const double H = gamma/(gamma-1.0)*pr/r + 0.5*(un*un+ut2);
            const double utg2 = max(0.0, 2.0*(H - gamma/(gamma-1.0)*pg/rg) - un*un);
            const double sc = ut2>1.0e-10 ? min(2.0,sqrt(utg2/ut2)) : 1.0;
            u = -f*un*c.nx + sc*tx;
            w = -f*un*c.nz + sc*tz;
            v = sc*v;
            r = rg;
            pr = pg;
        }
        else
        {
            u-=(1.0+f)*un*c.nx;
            v-=(1.0+f)*un*c.ny;
            w-=(1.0+f)*un*c.nz;
        }

        (*V[0])(c.i,c.j,c.k)=r;
        (*V[1])(c.i,c.j,c.k)=r*u;
        (*V[2])(c.i,c.j,c.k)=r*v;
        (*V[3])(c.i,c.j,c.k)=r*w;
        (*V[4])(c.i,c.j,c.k)=pr/(gamma-1.0)+0.5*r*(u*u+v*v+w*w);
    }

    // cells deep inside the bodies keep a harmless state
    for(size_t q=0; q<deep.size(); q+=3)
    for(int n=0; n<5; ++n)
    (*V[n])(deep[q],deep[q+1],deep[q+2]) = Qfs[n];
}

// MPI exchange, then physical boundaries: transmissive (zero gradient) or slip wall (mirror)
