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

#include"bc_ikepsilon_spalding.h"
#include"fdm.h"
#include"lexer.h"
#include"turbulence.h"

bc_ikepsilon_spalding::bc_ikepsilon_spalding(lexer* p)
    : wall_function_spalding(p)
{
}

bc_ikepsilon_spalding::~bc_ikepsilon_spalding()
{
}

void bc_ikepsilon_spalding::bckeps_start(fdm* a, lexer* p, field& kin, field& eps, int gcval, turbulence* pturb)
{
    int q;

    if(gcval==20 && p->B11!=0)
    {
        QGC4LOOP
        if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
        wall_nut(a, p, kin, eps, p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5], p->gcd4[q]);
        
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0)
            {
                a->rhsvec.V[n] -= a->M.s[n]*kin(i-1,j,k);
                a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
                a->rhsvec.V[n] -= a->M.n[n]*kin(i+1,j,k);
                a->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0 && p->j_dir==1)
            {
                a->rhsvec.V[n] -= a->M.e[n]*kin(i,j-1,k);
                a->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0 && p->j_dir==1)
            {
                a->rhsvec.V[n] -= a->M.w[n]*kin(i,j+1,k);
                a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
                a->rhsvec.V[n] -= a->M.b[n]*kin(i,j,k-1);
                a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
                a->rhsvec.V[n] -= a->M.t[n]*kin(i,j,k+1);
                a->M.t[n] = 0.0;
            }

            ++n;
        }
    }

    if(gcval==30 && p->B11!=0)
    {
        QGC4LOOP
        if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43 || (p->gcb4[q][4]==3 && p->gcb4[q][4]==6))
        wall_omega(a, p, kin, eps, p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5], p->gcd4[q]);
        
        n=0;
        LOOP
        {
            if(p->flag4[Im1JK]<0)
            {
                a->rhsvec.V[n] -= a->M.s[n]*eps(i-1,j,k);
                a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
                a->rhsvec.V[n] -= a->M.n[n]*eps(i+1,j,k);
                a->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0 && p->j_dir==1)
            {
                a->rhsvec.V[n] -= a->M.e[n]*eps(i,j-1,k);
                a->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0 && p->j_dir==1)
            {
                a->rhsvec.V[n] -= a->M.w[n]*eps(i,j+1,k);
                a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
                a->rhsvec.V[n] -= a->M.b[n]*eps(i,j,k-1);
                a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
                a->rhsvec.V[n] -= a->M.t[n]*eps(i,j,k+1);
                a->M.t[n] = 0.0;
            }

            ++n;
        }
    }
} 