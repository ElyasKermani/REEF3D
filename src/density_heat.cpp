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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"density_heat.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"heat_material.h"

density_heat::density_heat(lexer* p, heat *& ppheat)
{
    H=0.0;
    pheat = ppheat;

    if(p->j_dir==0)
    psi = p->F45*(1.0/2.0)*(p->DXM+p->DZM);

    if(p->j_dir==1)
    psi = p->F45*(1.0/3.0)*(p->DXM+p->DYM+p->DZM);
}

density_heat::~density_heat()
{
}

double density_heat::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double temp;
    double ro_1,ro_2;
	
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        
        temp = 0.5*(pheat->val(i,j,k) + pheat->val(i+aa,j+bb,k+cc));

        ro_1 = heat_material::density(p,1,temp);
        ro_2 = heat_material::density(p,2,temp);

        if(p->F30==0 && p->F80==0)
        H=1.0;

        else
        {
            if(phival>psi)
            H=1.0;

            if(phival<-psi)
            H=0.0;

            if(fabs(phival)<=psi)
            H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        }
        
        roval = ro_1*H + ro_2*(1.0-H);

	return roval;		
}
