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

#include"density_rheo.h"
#include"lexer.h"
#include"fdm.h"

density_rheo::density_rheo(lexer* p) 
{
}

density_rheo::~density_rheo()
{
}

double density_rheo::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{

    phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

    if(phival>p->psi)
    H=1.0;

    if(phival<-p->psi)
    H=0.0;

    if(fabs(phival)<=p->psi)
    H=0.5*(1.0 + phival/p->psi + (1.0/PI)*sin((PI*phival)/p->psi));
    
    roval = p->W1*H + p->W3*(1.0-H);

    return roval;          

}




