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

#include"lexer.h"
#include"fdm.h"
#include"sliceint.h"
#include"field.h"
#include"ghostcell.h"

void ghostcell::start1(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax(p,f,1);
	gcparacox(p,f,gcv);
    gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    
    // solid ghostcells
    starttime=timer();
	QQGC1LOOP
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
	gcdistro1(p,f,p->gcb1[qq][0], p->gcb1[qq][1], p->gcb1[qq][2], p->gcb1[qq][5], p->gcd1[qq], gcv, p->gcb1[qq][4], p->gcb1[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,1);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 1, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 1, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 1, 3);
    
    
    if(p->Y40==1  || p->Y40==3)
    dgcpol1(p,f,gcv);
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

void ghostcell::start2(lexer *p, field& f, int gcv)
{
    
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax(p,f,2);
	gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
    gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    
    if(p->j_dir==1)
    {
    starttime=timer();
	QQGC2LOOP
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
	gcdistro2(p,f,p->gcb2[qq][0], p->gcb2[qq][1], p->gcb2[qq][2], p->gcb2[qq][5], p->gcd2[qq], gcv, p->gcb2[qq][4], p->gcb2[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,2);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 2, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 2, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 2, 3);
    }
    
    if(p->Y40==1  || p->Y40==3)
    dgcpol2(p,f,gcv);
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

void ghostcell::start3(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax(p,f,3);
	gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
    gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    starttime=timer();
    QQGC3LOOP
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
	gcdistro3(p,f,p->gcb3[qq][0], p->gcb3[qq][1], p->gcb3[qq][2], p->gcb3[qq][5], p->gcd3[qq], gcv, p->gcb3[qq][4], p->gcb3[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,3);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 3, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 3, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 3, 3);
    
    
    if(p->Y40==1  || p->Y40==3)
    dgcpol3(p,f,gcv);
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

void ghostcell::start4(lexer *p, field &f, int gcv)
{
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax(p,f,4);
	gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
    gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
	
	starttime=timer();
	QQGC4LOOP
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
	gcdistro4(p,f,p->gcb4[qq][0],p->gcb4[qq][1], p->gcb4[qq][2], p->gcb4[qq][5], p->gcd4[qq], gcv, p->gcb4[qq][4], p->gcb4[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,4);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 4, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 4, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 4, 3);
    
    if(p->Y40==1  || p->Y40==3)
    dgcpol4(p,f,gcv);
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

void ghostcell::start4a(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax4a(p,f,5);
	gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    } 
    
    starttime=timer();
	QQGC4ALOOP
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
	gcdistro4a(p,f,p->gcb4a[qq][0], p->gcb4a[qq][1], p->gcb4a[qq][2], p->gcb4a[qq][5], p->gcd4a[qq], gcv, p->gcb4a[qq][4], p->gcb4a[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,4);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 4, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 4, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 4, 3);
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

void ghostcell::start4a_sum(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparax4a_sum(p,f,5);
	//gcparacox_sum(p,f,gcv);
	//gcparacox_sum(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    } 
    
    starttime=timer();
	QQGC4ALOOP
	gcdistro4a(p,f,p->gcb4a[qq][0], p->gcb4a[qq][1], p->gcb4a[qq][2], p->gcb4a[qq][5], p->gcd4a[qq], gcv, p->gcb4a[qq][4], p->gcb4a[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    // periodic ghostcells
    gcperiodicx(p,f,4);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 4, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 4, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 4, 3);
    
    
    if(p->M10>0)
	gcparacox(p,f,gcv);
}

