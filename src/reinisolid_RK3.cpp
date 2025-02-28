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

#include"reinisolid_RK3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf_rig.h"

reinisolid_RK3::reinisolid_RK3(lexer* p) : epsi(p->F45*p->DXM),f(p),frk1(p),frk2(p),L(p),dt(p)
{
	if(p->S50==1)
	gcval_topo=151;

	if(p->S50==2)
	gcval_topo=152;

	if(p->S50==3)
	gcval_topo=153;
	
	if(p->S50==4)
	gcval_topo=150;

	gcval_initopo=150;


	prdisc = new reinidisc_fsf_rig(p);

    time_preproc(p);    
}

reinisolid_RK3::~reinisolid_RK3()
{
}

void reinisolid_RK3::start(lexer *p, fdm *a, ghostcell *pgc, field &f)
{ 
    gcval=gcval_topo;

	pgc->start4a(p,f,gcval);
	
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<"initializing solid..."<<endl<<endl;
	reiniter=2*int(p->maxlength/(p->F43*p->DXM));
    gcval=gcval_initopo;
	pgc->start4a(p,f,gcval);
    }

	if(p->count>0)
	step(p,a);


    for(int q=0;q<reiniter;++q)
    {
	// Step 1
	prdisc->start(p,a,pgc,f,L,5);

	ALOOP
	frk1.V[IJK] = f.V[IJK] + dt.V[IJK]*L.V[IJK];

	pgc->start4a(p,frk1,gcval);
    
    
    // Step 2
    prdisc->start(p,a,pgc,frk1,L,5);

	ALOOP
	frk2.V[IJK]=  0.75*f.V[IJK] + 0.25*frk1.V[IJK] + 0.25*dt.V[IJK]*L.V[IJK];

	pgc->start4a(p,frk2,gcval);
    

    // Step 3
    prdisc->start(p,a,pgc,frk2,L,5);

	ALOOP
	f.V[IJK] = (1.0/3.0)*f.V[IJK] + (2.0/3.0)*frk2.V[IJK] + (2.0/3.0)*dt.V[IJK]*L.V[IJK];

	pgc->start4a(p,f,gcval);
	}
}

void reinisolid_RK3::step(lexer* p, fdm *a)
{
	reiniter=p->S37;
}

void reinisolid_RK3::time_preproc(lexer* p)
{	
    n=0;
	ALOOP
	{
	if(p->j_dir==0)
    dt.V[IJK] = p->F43*MIN(p->DXP[IP],p->DZP[KP]);
    
    if(p->j_dir==1)
	dt.V[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}
}


