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

#include"fluid_update_fsf_heat.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heat.h"
#include"heat_material.h"

fluid_update_fsf_heat::fluid_update_fsf_heat(lexer *p, fdm* a, ghostcell* pgc, heat *&ppheat)
{
    gcval_ro=1;
	gcval_visc=1;

	pheat = ppheat;
}

fluid_update_fsf_heat::~fluid_update_fsf_heat()
{
}

void fluid_update_fsf_heat::start(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
	double H=0.0;
	double temp;
    double visc_1,visc_2,ro_1,ro_2;
	p->volume1=0.0;
	p->volume2=0.0;

    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
    if(p->j_dir==0)        
    epsi = p->F45*(1.0/2.0)*(p->DXM+p->DZM);
        
    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DXM+p->DYM+p->DZM);

	LOOP
	{
        temp = pheat->val(i,j,k);

        ro_1   = heat_material::density(p,1,temp);
        ro_2   = heat_material::density(p,2,temp);
        visc_1 = heat_material::viscosity(p,1,temp);
        visc_2 = heat_material::viscosity(p,2,temp);

        if(p->F30==0 && p->F80==0)
        H=1.0;

        else
        {
            if(a->phi(i,j,k)>epsi)
            H=1.0;

            if(a->phi(i,j,k)<-epsi)
            H=0.0;

            if(fabs(a->phi(i,j,k))<=epsi)
            H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
        }

		a->ro(i,j,k)=      ro_1*H +   ro_2*(1.0-H);
		a->visc(i,j,k)= visc_1*H + visc_2*(1.0-H);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}

	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);

    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;

}

int fluid_update_fsf_heat::iocheck;
int fluid_update_fsf_heat::iter;
