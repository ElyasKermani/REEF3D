/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"bcmom_spalding.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

bcmom_spalding::bcmom_spalding(lexer* p)
    : surftens(p), wall_function_spalding(p), kappa(0.4)
{
    if(p->F50==1)
    gcval_phi=51;

    if(p->F50==2)
    gcval_phi=52;

    if(p->F50==3)
    gcval_phi=53;

    if(p->F50==4)
    gcval_phi=54;

    bckin=0;
    if(p->T10>0 || p->T10<20)
    bckin=1;
}

bcmom_spalding::~bcmom_spalding()
{
}

void bcmom_spalding::bcmom_start(fdm* a, lexer* p, ghostcell *pgc, turbulence *pturb, field& b, int gcval)
{
    int q;

    if(gcval==10 && p->B10!=0)
    {
        QGC1LOOP
        if((p->gcb1[q][4]==5 || p->gcb1[q][4]==21 || p->gcb1[q][4]==22 || p->gcb1[q][4]==41 || p->gcb1[q][4]==42 || p->gcb1[q][4]==43) && p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
        wall_tau_u(a, p, pturb, b, p->gcb1[q][0], p->gcb1[q][1], p->gcb1[q][2], p->gcb1[q][3], p->gcb1[q][4], p->gcd1[q]);
        
        QGCDF1LOOP
        wall_tau_u(a, p, pturb, b, p->gcdf1[q][0], p->gcdf1[q][1], p->gcdf1[q][2], p->gcdf1[q][3], p->gcdf1[q][4], 0.5*p->DXM);
    }

    if(gcval==11 && p->B10!=0 && p->j_dir==1)
    {
        QGC2LOOP
        if((p->gcb2[q][4]==5 || p->gcb2[q][4]==21 || p->gcb2[q][4]==22 || p->gcb2[q][4]==41 || p->gcb2[q][4]==42 || p->gcb2[q][4]==43) && p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
        wall_tau_v(a, p, pturb, b, p->gcb2[q][0], p->gcb2[q][1], p->gcb2[q][2], p->gcb2[q][3], p->gcb2[q][4], p->gcd2[q]);
        
        QGCDF2LOOP
        wall_tau_v(a, p, pturb, b, p->gcdf2[q][0], p->gcdf2[q][1], p->gcdf2[q][2], p->gcdf2[q][3], p->gcdf2[q][4], 0.5*p->DXM);
    }

    if(gcval==12 && p->B10!=0)
    {
        QGC3LOOP
        if((p->gcb3[q][4]==5 || p->gcb3[q][4]==21 || p->gcb3[q][4]==22 || p->gcb3[q][4]==41 || p->gcb3[q][4]==42 || p->gcb3[q][4]==43) && p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
        wall_tau_w(a, p, pturb, b, p->gcb3[q][0], p->gcb3[q][1], p->gcb3[q][2], p->gcb3[q][3], p->gcb3[q][4], p->gcd3[q]);
        
        QGCDF3LOOP
        wall_tau_w(a, p, pturb, b, p->gcdf3[q][0], p->gcdf3[q][1], p->gcdf3[q][2], p->gcdf3[q][3], p->gcdf3[q][4], 0.5*p->DXM);
    }
    
    // Apply surface tension effects same as in original bcmom
    surface_tension(a, p, a->phi, gcval);
} 