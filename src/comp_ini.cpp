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

void comp_euler::read_input()
{
    ifstream in("comp.txt");

    if(!in)
    {
        if(p->mpirank==0)
        cout<<"REEF3D::COMP: comp.txt not found, using the Sod shock tube in x"<<endl;
        return;
    }

    string line;
    while(getline(in,line))
    {
        istringstream ss(line);
        string key;
        if(!(ss>>key) || key[0]=='#')
        continue;

        if(key=="gamma") ss>>gamma;
        else if(key=="cfl") { ss>>cfl; p->N47=cfl; }
        else if(key=="tend") { double t; ss>>t; p->N41=t; }
        else if(key=="recon") ss>>recon;
        else if(key=="flux") ss>>flux;
        else if(key=="dir") ss>>dir;
        else if(key=="x0") ss>>x0;
        else if(key=="left") ss>>rhoL>>uL>>pL;
        else if(key=="right") ss>>rhoR>>uR>>pR;
        else if(key=="case")
        {
            string c;
            ss>>c;
            casetype = (c=="shuosher") ? 1 : (c=="freestream" ? 2 : (c=="shearwave" ? 3 : (c=="channelwave" ? 4 : (c=="entropywave" ? 5 : (c=="blast" ? 6 : 0)))));
        }
        else if(key=="bc") ss>>bc[0]>>bc[1]>>bc[2]>>bc[3]>>bc[4]>>bc[5];
        else if(key=="zwall_from") ss>>zwall_from;
        else if(key=="curvature") ss>>curvature;
        else if(key=="weno_nonuniform") ss>>nuweno;
        else if(key=="viscosity") ss>>mu;
        else if(key=="shear") ss>>shearU>>shearH;
        else if(key=="channel") ss>>chanU>>chanH>>chanAng>>chanN0;
        else if(key=="blast") ss>>blastE>>blastR>>blastX>>blastZ>>blastP0>>blastRho0;
        else if(key=="sponge") ss>>sponge_in>>sponge_out>>sponge_sigma;
        else if(key=="sponge_mode") ss>>sponge_mode;
        else if(key=="normal")
        {
            ss>>nrm[0]>>nrm[1]>>nrm[2];
            usenormal=1;
            p->A830=1;
            p->A830_nx=nrm[0]; p->A830_ny=nrm[1]; p->A830_nz=nrm[2];
        }
        else if(key=="restart") ss>>restart;
        else if(key=="sutherland") ss>>sutherland;
        else if(key=="corner_fix") ss>>corner_x>>corner_z>>corner_n;
        else if(key=="top_state")
        {
            // top_state x0 rho V angle_deg p: the +z boundary (type 2) takes this state for x >= x0
            ss>>top_x0>>p->A828_rho>>p->A828_V>>p->A828_ang>>p->A828_p;
            p->A828_x0=top_x0;
        }
        else if(key=="prandtl") ss>>prandtl;
        else if(key=="noslip") ss>>noslip;
        else if(key=="solid_poly")
        {
            // solid_poly n x1 z1 ... xn zn: a polygon in the x-z plane (any orientation)
            int nv; ss>>nv; std::vector<double> v(2*nv);
            for(int q=0; q<2*nv; ++q) ss>>v[q];
            polys.push_back(v);
        }
        else if(key=="solid_box")
        {
            double a,b,c,d; ss>>a>>b>>c>>d;
            polys.push_back({a,c, b,c, b,d, a,d});
        }
        else if(key=="freestream")
        {
            // rho, velocity magnitude, angle in the x-z plane (degrees, positive upwards), pressure
            ss>>p->A813_rho>>p->A813_V>>p->A813_ang>>p->A813_p;
        }
        else if(key=="print")
        {
            double t;
            while(ss>>t)
            prints.push_back(t);
        }
    }
    sort(prints.begin(),prints.end());
}

void comp_euler::riemann_state(double s, double sl, double, double *Q) const
{
    // s: coordinate along the problem direction; sl unused here
    double rho,un,pr;

    if(casetype==2)
    {
        for(int n=0; n<5; ++n)
        Q[n]=Qfs[n];
        return;
    }

    if(casetype==1)   // Shu-Osher: Mach 3 shock into a sinusoidal density field, x0 = -4
    {
        if(s<x0)
        {
            rho=3.857143; un=2.629369; pr=10.33333;
        }
        else
        {
            rho=1.0+0.2*sin(5.0*s); un=0.0; pr=1.0;
        }
    }
    else
    {
        if(s<x0)
        {
            rho=rhoL; un=uL; pr=pL;
        }
        else
        {
            rho=rhoR; un=uR; pr=pR;
        }
    }

    double vel[3]={0.0,0.0,0.0};
    if(usenormal)
    for(int q=0; q<3; ++q) vel[q]=un*nrm[q];
    else
    vel[dir]=un;

    Q[0]=rho;
    Q[1]=rho*vel[0];
    Q[2]=rho*vel[1];
    Q[3]=rho*vel[2];
    Q[4]=pr/(gamma-1.0) + 0.5*rho*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
}

void comp_euler::initialize()
{
    double Q[5];
    std::vector<double> Qa(5);

    if(casetype==6 && restart.empty())
    {
        // point blast: ambient state, and the energy E spread uniformly over the fluid cells within r0
        double vol=0.0;
    ILOOP
    JLOOP
    KLOOP
        {
            const double dxr=p->XP[IP]-blastX, dzr=p->ZP[KP]-blastZ;
            if(sqrt(dxr*dxr+dzr*dzr)<blastR && !solid(i,j,k))
            vol+=p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        }
        vol=pgc->globalsum(vol);
    ILOOP
    JLOOP
    KLOOP
        {
            const double dxr=p->XP[IP]-blastX, dzr=p->ZP[KP]-blastZ;
            const bool src = sqrt(dxr*dxr+dzr*dzr)<blastR && !solid(i,j,k) && vol>0.0;
            // E is per unit length in y for a 2D run: divide by the slab thickness
            const double ly = p->gknoy>1 ? 1.0 : p->DYN[JP];
            const double eint = blastP0/(gamma-1.0) + (src ? blastE*ly/vol : 0.0);
            (*U[0])(i,j,k)=blastRho0; (*U[1])(i,j,k)=0.0; (*U[2])(i,j,k)=0.0; (*U[3])(i,j,k)=0.0;
            (*U[4])(i,j,k)=eint;
        }
        if(p->mpirank==0)
        cout<<"REEF3D::COMP blast: energy "<<blastE<<" in "<<vol<<" (volume of the source cells)"<<endl;
        return;
    }

    if(!restart.empty())
    {
        ifstream in(restart.c_str());
        string line;
        getline(in,line);
        { istringstream h(line); string a,b; h>>a>>b>>restart_time; }
        getline(in,line);
    ILOOP
    JLOOP
    KLOOP
        {
            double x,y,z,r,u,v,w,pr;
            in>>x>>y>>z>>r>>u>>v>>w>>pr;
            (*U[0])(i,j,k)=r; (*U[1])(i,j,k)=r*u; (*U[2])(i,j,k)=r*v; (*U[3])(i,j,k)=r*w;
            (*U[4])(i,j,k)=pr/(gamma-1.0)+0.5*r*(u*u+v*v+w*w);
        }
        p->simtime=restart_time;
        if(p->mpirank==0)
        cout<<"REEF3D::COMP restart from "<<restart<<" at t = "<<restart_time<<endl;
        return;
    }

    // initial cell values are cell averages (3-point Gauss rule per direction) of the point states,
    // which keeps smooth initial data at the scheme's order and places a discontinuity inside its cell
    const double gq[3]={-0.7745966692414834,0.0,0.7745966692414834}, gw[3]={5.0/18.0,8.0/18.0,5.0/18.0};
    auto point=[&](double x, double y, double z, double *Qp)
    {
        double s = dir==0 ? x : (dir==1 ? y : z);
        if(usenormal) s = nrm[0]*x + nrm[1]*y + nrm[2]*z;
        riemann_state(s,0.0,0.0,Qp);
        if(casetype==4)
        {
            const double a=chanAng*PI/180.0;
            const double nn=-x*sin(a)+z*cos(a)-chanN0;
            const double u = (nn>0.0 && nn<chanH) ? chanU*sin(PI*nn/chanH) : 0.0;
            const double r=1.0, pr=1.0/gamma;
            Qp[0]=r; Qp[1]=r*u*cos(a); Qp[2]=0.0; Qp[3]=r*u*sin(a); Qp[4]=pr/(gamma-1.0)+0.5*r*u*u;
        }
        if(casetype==3)
        {
            const double u=shearU*cos(PI*z/shearH), r=1.0, pr=1.0/gamma;
            Qp[0]=r; Qp[1]=r*u; Qp[2]=0.0; Qp[3]=0.0; Qp[4]=pr/(gamma-1.0)+0.5*r*u*u;
        }
    };

    ILOOP
    JLOOP
    KLOOP
    {
        double s = dir==0 ? p->XP[IP] : (dir==1 ? p->YP[JP] : p->ZP[KP]);
        if(usenormal)
        s = nrm[0]*p->XP[IP] + nrm[1]*p->YP[JP] + nrm[2]*p->ZP[KP];
        riemann_state(s,0.0,0.0,Q);

        if(casetype==0 || casetype==3 || casetype==4)
        {
            const double xc=p->XP[IP], yc=p->YP[JP], zc=p->ZP[KP];
            const double hx=0.5*p->DXN[IP], hy=0.5*p->DYN[JP], hz=0.5*p->DZN[KP];
            const int nx = p->gknox>1 ? 3 : 1, ny = p->gknoy>1 ? 3 : 1, nz = p->gknoz>1 ? 3 : 1;
            for(int n=0; n<5; ++n) Qa[n]=0.0;
            for(int a=0; a<nx; ++a)
            for(int b=0; b<ny; ++b)
            for(int c=0; c<nz; ++c)
            {
                const double w=(nx>1?gw[a]:1.0)*(ny>1?gw[b]:1.0)*(nz>1?gw[c]:1.0);
                double Qp[5];
                point(xc+(nx>1?gq[a]*hx:0.0), yc+(ny>1?gq[b]*hy:0.0), zc+(nz>1?gq[c]*hz:0.0), Qp);
                for(int n=0; n<5; ++n) Qa[n]+=w*Qp[n];
            }
            for(int n=0; n<5; ++n) Q[n]=Qa[n];
            for(int n=0; n<5; ++n) (*U[n])(i,j,k)=Q[n];
            continue;
        }

        if(casetype==4)
        {
            // channel at angle a (degrees): along t = (cos a, sin a), across n = (-sin a, cos a); between the walls
            // n0 < n < n0 + H the velocity is U sin(pi (n - n0)/H) t, decaying exactly as exp(-nu (pi/H)^2 t)
            const double a=chanAng*PI/180.0;
            const double nn=-p->XP[IP]*sin(a)+p->ZP[KP]*cos(a)-chanN0;
            const double u = (nn>0.0 && nn<chanH) ? chanU*sin(PI*nn/chanH) : 0.0;
            const double r=1.0, pr=1.0/gamma;
            Q[0]=r; Q[1]=r*u*cos(a); Q[2]=0.0; Q[3]=r*u*sin(a); Q[4]=pr/(gamma-1.0)+0.5*r*u*u;
        }

        if(casetype==5)   // entropy wave rho = 1 + 0.2 sin(2 pi x), u = 1, p = 1/gamma: exact solution translates
        {
            // cell average of the density (the scheme works with cell averages)
            const double xl=p->XN[IP], xr=p->XN[IP1], k2=2.0*PI;
            const double r=1.0+0.2*(cos(k2*xl)-cos(k2*xr))/(k2*(xr-xl)), u=1.0, pr=1.0/gamma;
            Q[0]=r; Q[1]=r*u; Q[2]=0.0; Q[3]=0.0; Q[4]=pr/(gamma-1.0)+0.5*r*u*u;
        }

        if(casetype==3)   // decaying shear wave u = U cos(pi z/H): exact decay exp(-nu (pi/H)^2 t)
        {
            const double u=shearU*cos(PI*p->ZP[KP]/shearH), r=1.0, pr=1.0/gamma;
            Q[0]=r; Q[1]=r*u; Q[2]=0.0; Q[3]=0.0; Q[4]=pr/(gamma-1.0)+0.5*r*u*u;
        }

        for(int n=0; n<5; ++n)
        (*U[n])(i,j,k)=Q[n];
    }
}

