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

#include"wave_lib_wcp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_wcp::wave_lib_wcp(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    // read header
    read_header(p,pgc);
    allocate(p,pgc);
    
    
    // time_interpol
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: wave coupling FNPF->CFD "<<endl;
    cout<<p->mpirank<<" WCP Nx: "<<Nx<<" Ny: "<<Ny<<" Nz: "<<Nz<<" . jdir: "<<jdir<<endl;
    cout<<p->mpirank<<" WCP Xs: "<<Xstart<<" Xe: "<<Xend<<" Ys: "<<Ystart<<" Ye: "<<Yend<<endl;
    cout<<p->mpirank<<" WCP numiter: "<<numiter<<" t_start: "<<t_start<<" t_end: "<<t_end<<endl;
    }
        
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
    
    startup=0;
    endseries=0;
}

wave_lib_wcp::~wave_lib_wcp()
{
}

double wave_lib_wcp::wave_u(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(endseries==0)
    vel = space_interpol(p,U,x,y,z);

    return cosgamma*vel;
}

double wave_lib_wcp::wave_v(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(endseries==0)
    vel = space_interpol(p,V,x,y,z);

    return singamma*vel;
}

double wave_lib_wcp::wave_w(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(endseries==0)
    vel = space_interpol(p,W,x,y,z);

    return vel;
}

double wave_lib_wcp::wave_eta(lexer *p, double x, double y)
{
    double eta=0.0;
    
    if(endseries==0)
    eta = plane_interpol(p,E,x,y);

    return eta;
}

double wave_lib_wcp::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    return fi;
}

void wave_lib_wcp::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_wcp::wave_prestep(lexer *p, ghostcell *pgc)
{
    // only at startup
    if(startup==0)
    {
        t1 = simtime[0];
        t2 = simtime[1];
        q1 = 0;
        q2 = 1;
    
        read_result(p,pgc,E1,U1,V1,W1,q1);
        read_result(p,pgc,E2,U2,V2,W2,q2);
        startup=1;
    }
    
    q1n=q1;
    q2n=q2;
    
    // check: open next timestep
           
    // find q1
    while(simtime[q1+1]<=p->simtime+p->I241)
    ++q1;
        
    // find q2
    while(simtime[q2]<p->simtime+p->I241)
    ++q2;
    
    if(q2>=numiter)
    endseries=1;
        
    q1=MIN(q1,numiter);
    q2=MIN(q2,numiter);
    
    if(q1==q2)
    ++q2;
    
    
    
        if(q1!=q1n)
        {
        // Open File 1
        filename(p,pgc,q1);
        read_result(p,pgc,E1,U1,V1,W1,q1);
        }
        
        
        if(q2!=q2n)
        {
        // Open File 2
        filename(p,pgc,q2);
        read_result(p,pgc,E2,U2,V2,W2,q2);
        }
        

        deltaT = simtime[q2]-simtime[q1];
        
        if(p->mpirank==0)
        cout<<"WCP  q1: "<<q1<<" q2: "<<q2<<" deltaT: "<<deltaT<<" simtime[q1]: "<<simtime[q1]<<" simtime[q2]: "<<simtime[q2]<<endl;

        t1 = (simtime[q2]-(p->simtime+p->I241))/deltaT;
        t2 = ((p->simtime+p->I241)-simtime[q1])/deltaT;
        
        
    // time interpolation
    if(q1!=q1n || q2!=q2n)
    time_interpol(p);
    
    
}