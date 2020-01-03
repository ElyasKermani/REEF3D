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

#include"idiff2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

idiff2::idiff2(lexer* p):gradient(p)
{
	vfm=vft=0.0;
	
	if(p->D22==1)
	vfm=1.0;
	
	if(p->D23==1)
	vft=1.0;
}

idiff2::~idiff2()
{
}

void idiff2::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
     count=0;
	 
	 double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;

     sqd = (1.0/(p->dx*p->dx));
     
    ULOOP
	{
	ev_ijk=a->eddyv(i,j,k);
	ev_ip_j_k=a->eddyv(i+1,j,k);
	ev_i_jm_k=a->eddyv(i,j-1,k);
	ev_i_jp_k=a->eddyv(i,j+1,k);
	ev_i_j_km=a->eddyv(i,j,k-1);
	ev_i_j_kp=a->eddyv(i,j,k+1);
	
	visc_ijk=a->visc(i,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_jm_k=a->visc(i,j-1,k);
	visc_i_jp_k=a->visc(i,j+1,k);
	visc_i_j_km=a->visc(i,j,k-1);
	visc_i_j_kp=a->visc(i,j,k+1);
	
	visc_ddy_p = (vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k + vfm*visc_i_jp_k+visc_i_jp_k + vfm*a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddy_m = (vfm*visc_i_jm_k+ev_i_jm_k  +vfm*a->visc(i+1,j-1,k)+a->eddyv(i+1,j-1,k) + vfm*visc_ijk+ev_ijk + vfm*a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	visc_ddz_p = (vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k + vfm*visc_i_j_kp+ev_i_j_kp + vfm*a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddz_m = (vfm*a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + vfm*a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + vfm*visc_ijk+ev_ijk + vfm*visc_ip_j_k+ev_ip_j_k)*0.25;
		
        
	a->M.p[count] +=  2.0*(vfm*visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP]*p->DXP[IP])
				   + 2.0*(vfm*visc_ijk+ev_ijk)/(p->DXN[IM1]*p->DXP[IP])
				   + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				   + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				   + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				   + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				   + CPOR1/(alpha*p->dt);
				  
	a->rhsvec.V[count] += ((v(i+1,j,k)-v(i,j,k))*visc_ddy_p - (v(i+1,j-1,k)-v(i,j-1,k))*visc_ddy_m)/(p->DXP[IP]*p->DYN[JP])
						 + ((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP])

						 + (CPOR1*u(i,j,k))/(alpha*p->dt);
	 
	 a->M.s[count] -= 2.0*(vfm*visc_ijk+ev_ijk)/(p->DXN[IM1]*p->DXP[IP]);
	 a->M.n[count] -= 2.0*(vfm*visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP]*p->DXP[IP]);
	 
	 a->M.e[count] -= visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] -= visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] -= visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] -= visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
	 
	 ++count;
	}
}


