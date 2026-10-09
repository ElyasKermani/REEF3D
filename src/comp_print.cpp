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

void comp_euler::print_lagoon(int count)
{
    if(!plagoon)
    plagoon = new lagoon_output(p,pgc,"COMP");

    const int nx=p->knox+1, ny=p->knoy+1, nz=p->knoz+1;
    const size_t n=size_t(nx)*ny*nz;
    p->pointnum=int(n);

    std::vector<float> pts(3*n), vel(3*n), rho(n), pre(n), mach(n), lsv(n);
    size_t q=0;
    for(int kk=0; kk<nz; ++kk)
    for(int jj=0; jj<ny; ++jj)
    for(int ii=0; ii<nx; ++ii)
    {
        double s[5]={0,0,0,0,0}, w=0.0, ls=0.0;
        for(int a=ii-1; a<=ii; ++a)
        for(int b=jj-1; b<=jj; ++b)
        for(int c=kk-1; c<=kk; ++c)
        {
            // at a physical boundary only the cells inside the domain count (corner ghosts are not filled)
            if((a<0 && p->origin_i==0) || (a>=p->knox && p->origin_i+p->knox==p->gknox)) continue;
            if((b<0 && p->origin_j==0) || (b>=p->knoy && p->origin_j+p->knoy==p->gknoy)) continue;
            if((c<0 && p->origin_k==0) || (c>=p->knoz && p->origin_k+p->knoz==p->gknoz)) continue;
            for(int m=0; m<5; ++m)
            s[m]+=(*U[m])(a,b,c);
            ls+=(*phi)(a,b,c);
            w+=1.0;
        }
        for(int m=0; m<5; ++m)
        s[m]/=w;

        double r,u,v,ww,pr;
        prim(s,r,u,v,ww,pr);
        pts[3*q]=float(p->XN[ii+marge]); pts[3*q+1]=float(p->YN[jj+marge]); pts[3*q+2]=float(p->ZN[kk+marge]);
        vel[3*q]=float(u); vel[3*q+1]=float(v); vel[3*q+2]=float(ww);
        rho[q]=float(r); pre[q]=float(pr); lsv[q]=float(ls/w);
        mach[q]=float(sqrt(u*u+v*v+ww*ww)/sqrt(gamma*max(pr,1.0e-14)/r));
        ++q;
    }

    std::ostringstream h;
    const char *names[5]={"velocity","pressure","density","Mach","levelset"};
    const int comps[5]={3,1,1,1,1};
    long long off=0;
    h<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
    h<<"<Piece NumberOfPoints=\""<<n<<"\" NumberOfCells=\"0\">\n<Points>\n";
    h<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<off<<"\"/>\n</Points>\n<PointData>\n";
    off+=sizeof(int)+3*n*sizeof(float);
    for(int a=0; a<5; ++a)
    {
        h<<"<DataArray type=\"Float32\" Name=\""<<names[a]<<"\"";
        if(comps[a]>1) h<<" NumberOfComponents=\""<<comps[a]<<"\"";
        h<<" format=\"appended\" offset=\""<<off<<"\"/>\n";
        off+=sizeof(int)+comps[a]*n*sizeof(float);
    }
    h<<"</PointData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";

    const std::string head=h.str();
    std::vector<char> buf(head.begin(),head.end());
    const size_t data_start=buf.size();
    auto put=[&](const std::vector<float> &v)
    {
        int bytes=int(v.size()*sizeof(float));
        const char *c=reinterpret_cast<const char*>(&bytes);
        buf.insert(buf.end(),c,c+sizeof(int));
        const char *d=reinterpret_cast<const char*>(v.data());
        buf.insert(buf.end(),d,d+bytes);
    };
    put(pts); put(vel); put(pre); put(rho); put(mach); put(lsv);

    plagoon->vtu_piece(p,pgc,buf,data_start,count);
}

void comp_euler::print(int count)
{
    if(p->P18>0)
    {
        boundary(U);
        print_lagoon(count);
        if(p->P18==1)
        return;
    }

    if(p->mpirank==0)
    mkdir("./REEF3D_COMP",0777);

    pgc->globalmax(0.0);   // barrier

    char name[200];
    snprintf(name,sizeof(name),"./REEF3D_COMP/comp-%04i-%05i.dat",count,p->mpirank);
    ofstream out(name);
    out<<setprecision(12);
    out<<"# time "<<p->simtime<<" steps "<<p->count<<endl;
    out<<"# x y z rho u v w p"<<endl;

    ILOOP
    JLOOP
    KLOOP
    {
        double Q[5], r,u,v,w,pr;
        for(int n=0; n<5; ++n)
        Q[n]=(*U[n])(i,j,k);
        prim(Q,r,u,v,w,pr);

        out<<p->XP[IP]<<" "<<p->YP[JP]<<" "<<p->ZP[KP]<<" "<<r<<" "<<u<<" "<<v<<" "<<w<<" "<<pr<<endl;
    }
}

