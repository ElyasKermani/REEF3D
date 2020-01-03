/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"bc_ikepsilon.h"
#include"fdm.h"
#include"lexer.h"

bc_ikepsilon::bc_ikepsilon(lexer* p):roughness(p),kappa(0.4)
{
	deltax=p->dx;
}

bc_ikepsilon::~bc_ikepsilon()
{
}

void bc_ikepsilon::bckeps_start(fdm* a,lexer* p,field& kin,field& eps,int gcval)
{
	int q;

	if(gcval==20)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21  || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
		wall_law_kin(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
	}

	if(gcval==30)
	{
		QGC4LOOP
		if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
		wall_law_eps(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
	}


}
// ****************************
// WALL KIN
// ****************************
void bc_ikepsilon::wall_law_kin(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    double uvel,vvel,wvel;
    dist=0.5*p->DXM;

	i=ii;
	j=jj;
	k=kk;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

        pip=1;
        uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
        pip=0;

        pip=2;
        vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
        pip=0;

        pip=3;
        wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
        pip=0;

        u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);

		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));

	tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);


	a->M.p[id] += (pow(p->cmu,0.75)*pow(fabs(kin(i,j,k)),0.5)*uplus)/dist;
	a->rhsvec.V[id] += (tau*u_abs)/dist;

}

void bc_ikepsilon::wall_law_eps(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
	i=ii;
	j=jj;
	k=kk;
    dist=0.5*p->DXM;

	eps_star = (pow(p->cmu, 0.75)*pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),1.5)) / (0.4*dist);

	eps(i,j,k) = eps_star;
}



