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

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void patchBC::patchBC_gcb_convert(lexer *p, ghostcell *pgc)
{
    // convert gcbs
    int istart,iend,jstart,jend,kstart,kend,qn;
	
    int count=0;
    for(qn=0;qn<p->B221;++qn)
    {
        istart = p->posc_i(p->B221_xs[qn]);
        iend = p->posc_i(p->B221_xe[qn]);
        
        jstart = p->posc_j(p->B221_ys[qn]);
        jend = p->posc_j(p->B221_ye[qn]);
        
        kstart = p->posc_k(p->B221_zs[qn]);
        kend = p->posc_k(p->B221_ze[qn]);
        
        
        for(n=0;n<p->gcb4_count;++n)
		{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
			if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && p->gcb4[n][3]==p->B221_face[qn] && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
			{
			++count;
			p->gcb4[n][4]=21;
			}
		}
    }
    
    
    for(qn=0;qn<p->B231;++qn)
    {
        istart = p->posc_i(p->B231_xs[qn]);
        iend = p->posc_i(p->B231_xe[qn]);
        
        jstart = p->posc_j(p->B231_ys[qn]);
        jend = p->posc_j(p->B231_ye[qn]);
        
        kstart = p->posc_k(p->B231_zs[qn]);
        kend = p->posc_k(p->B231_ze[qn]);
        
        
        for(n=0;n<p->gcb4_count;++n)
		{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
			if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && p->gcb4[n][3]==p->B231_face[qn]  && (p->gcb4[n][4]==21||p->gcb4[n][4]==22))
			{
			++count;
			p->gcb4[n][4]=31;
			}
		}
    }
    
    
    // gcin / gcout


} 

