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

double comp_euler::timestep()
{
    double dtmin=1.0e20;
    const int gk[3]={p->gknox,p->gknoy,p->gknoz};

    ILOOP
    JLOOP
    KLOOP
    {
        double Q[5], r,u,v,w,pr;
        for(int n=0; n<5; ++n)
        Q[n]=(*U[n])(i,j,k);

        if(solid(i,j,k))
        continue;

        prim(Q,r,u,v,w,pr);
        const double c=sqrt(gamma*max(pr,1.0e-14)/r);

        // convective and viscous limits together, summed over the directions (stable for the
        // dimension-by-dimension update in 2D and 3D): 1/dt = sum (|u|+c)/h + 2 nu_eff sum 1/h^2
        const double nu = mu>0.0 ? muT(pr/r)/r*max(4.0/3.0,gamma/prandtl) : 0.0;
        double inv=0.0;
        if(gk[0]>1) inv+=(fabs(u)+c)/p->DXN[IP] + 2.0*nu/(p->DXN[IP]*p->DXN[IP]);
        if(gk[1]>1) inv+=(fabs(v)+c)/p->DYN[JP] + 2.0*nu/(p->DYN[JP]*p->DYN[JP]);
        if(gk[2]>1) inv+=(fabs(w)+c)/p->DZN[KP] + 2.0*nu/(p->DZN[KP]*p->DZN[KP]);
        if(inv>0.0) dtmin=min(dtmin,1.0/inv);
    }

    return p->N47*pgc->globalmin(dtmin);
}

int comp_euler::check_stage(field4 **V)
{
    int nbad=0;
    const int gk[3]={p->gknox,p->gknoy,p->gknoz};

    ILOOP
    JLOOP
    KLOOP
    {
        double Q[5], r,u,v,w,pr;
        for(int n=0; n<5; ++n)
        Q[n]=(*V[n])(i,j,k);
        prim(Q,r,u,v,w,pr);

        if(solid(i,j,k))
        continue;

        if(!(r>0.0) || !(pr>0.0))
        {
            ++nbad;
            (*B)(i,j,k)=1.0;
            if(gk[0]>1) { (*B)(i-1,j,k)=1.0; (*B)(i+1,j,k)=1.0; }
            if(gk[1]>1) { (*B)(i,j-1,k)=1.0; (*B)(i,j+1,k)=1.0; }
            if(gk[2]>1) { (*B)(i,j,k-1)=1.0; (*B)(i,j,k+1)=1.0; }
        }
    }
    return int(pgc->globalsum(double(nbad)));
}

