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

#include"density_heat_Bouss.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"

density_heat_Bouss::density_heat_Bouss(lexer* p, heat *&ppheat)
{
    if(p->j_dir==0)
    epsi = p->F45*(1.0/2.0)*(p->DXM+p->DZM);

    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DXM+p->DYM+p->DZM);

    if(p->H9==1)
    {
    T0_1 = p->H50_1 + 273.0;
    T0_2 = p->H50_2 + 273.0;
    }

    if(p->H9==2)
    {
    T0_1 = p->H50_2 + 273.0;
    T0_2 = p->H50_1 + 273.0;
    }

	pheat = ppheat;
}

density_heat_Bouss::~density_heat_Bouss()
{
}

double density_heat_Bouss::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double temp;
    double phival_i = a->phi(i,j,k);
    double phival_ip = a->phi(i+aa,j+bb,k+cc);
    double temp_i = pheat->val(i,j,k) + 273.0;
    double temp_ip = pheat->val(i+aa,j+bb,k+cc) + 273.0;

    phival = 0.5*(phival_i + phival_ip);
    temp = 0.5*(temp_i + temp_ip);

    if(p->H4==0)
    {
        if(p->H9==1)
        {
        ro_1 = p->W1 - p->W1*(temp - T0_1)/T0_1;
        ro_2 = p->W3 - p->W3*(temp - T0_2)/T0_2;
        }

        if(p->H9==2)
        {
        ro_1 = p->W3 - p->W3*(temp - T0_2)/T0_2;
        ro_2 = p->W1 - p->W1*(temp - T0_1)/T0_1;
        }
    }

    if(p->H4==1)
    {
        if(p->H9==1)
        {
        ro_1 = p->W1 - p->W1*p->H4_beta1*(temp - T0_1);
        ro_2 = p->W3 - p->W3*p->H4_beta2*(temp - T0_2);
        }

        if(p->H9==2)
        {
        ro_1 = p->W3 - p->W3*p->H4_beta2*(temp - T0_2);
        ro_2 = p->W1 - p->W1*p->H4_beta1*(temp - T0_1);
        }
    }

    if(phival>epsi)
    H=1.0;

    if(phival<-epsi)
    H=0.0;

    if(fabs(phival)<=epsi)
    H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));

    roval = ro_1*H + ro_2*(1.0-H);

	return roval;
}
