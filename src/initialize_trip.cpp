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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<iostream>

void initialize::trip_ini(lexer* p, fdm *a, ghostcell* pgc)
{
    const double A = p->I59;
    const double lx = p->xcoormax - p->xcoormin;
    const double ly = p->ycoormax - p->ycoormin;
    const double h = (p->F60>-1.0e20) ? p->F60 : (p->zcoormax - p->zcoormin);

    if(p->mpirank==0)
    cout<<"initial 3D trip  I 59 "<<A<<"  water depth "<<h<<endl;

    ULOOP
    {
        const double phi = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        const double z = p->ZP[KP] - p->zcoormin;
        if(phi>=0.0 && z>0.0 && z<h)
        {
            const double y = p->YP[JP] - p->ycoormin;
            const double fz = sin(PI*z/h);
            a->u(i,j,k) += 0.5*A*sin(2.0*PI*y/ly)*fz;
        }
    }

    VLOOP
    {
        const double phi = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
        const double z = p->ZP[KP] - p->zcoormin;
        if(phi>=0.0 && z>0.0 && z<h)
        {
            const double x = p->XP[IP] - p->xcoormin;
            const double y = p->YP[JP] - p->ycoormin;
            const double fz = sin(PI*z/h);
            a->v(i,j,k) += A*sin(2.0*PI*x/lx)*cos(2.0*PI*y/ly)*fz;
            a->v(i,j,k) += 0.5*A*sin(4.0*PI*x/lx)*sin(2.0*PI*y/ly)*sin(2.0*PI*z/h);
        }
    }

    WLOOP
    {
        const double phi = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        const double z = p->ZP[KP] - p->zcoormin;
        if(phi>=0.0 && z>0.0 && z<h)
        {
            const double x = p->XP[IP] - p->xcoormin;
            const double y = p->YP[JP] - p->ycoormin;
            const double fz = sin(PI*z/h);
            a->w(i,j,k) += A*cos(2.0*PI*x/lx)*sin(2.0*PI*y/ly)*fz;
        }
    }

    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
    pgc->start3(p,a->w,12);
}
