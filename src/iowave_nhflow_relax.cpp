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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::WL_relax(lexer *p, ghostcell *pgc, slice &WL, slice &depth)
{
    starttime=pgc->timer();
    
	count=0;
    SLICELOOP4
    {
		dg = distgen(p);
		db = distbeach(p);
        

		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            WETDRYDEEP
            WL(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*(eta(i,j) + depth(i,j)) + relax4_wg(i,j) * WL(i,j);
            ++count;
            }
		}
        
		
		// Numerical Beach
		if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            WL(i,j) = (1.0-relax4_nb(i,j))*depth(i,j) + relax4_nb(i,j)*WL(i,j);
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::U_relax(lexer *p, ghostcell *pgc, double *U, double *UH)
{
    starttime=pgc->timer();
    
    count=0;
    LOOP
    {
         dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            U[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*uval[count] + relax4_wg(i,j)*U[IJK];
            UH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*UHval[count] + relax4_wg(i,j)*UH[IJK];
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            U[IJK] = relax4_nb(i,j)*U[IJK];
            UH[IJK] = relax4_nb(i,j)*UH[IJK];
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::V_relax(lexer *p, ghostcell *pgc, double *V, double *VH)
{ 
    starttime=pgc->timer();
    
    count=0;
    if(p->j_dir==1)
    LOOP
    {
        dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            V[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*vval[count] + relax4_wg(i,j)*V[IJK];
            VH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*VHval[count] + relax4_wg(i,j)*VH[IJK];
            }
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1||p->B99==2||beach_relax==1)
		{	
            // Zone 2
            if(db<1.0e20)
            {
            V[IJK] = relax4_nb(i,j)*V[IJK];
            VH[IJK] = relax4_nb(i,j)*VH[IJK];
            }
            
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::W_relax(lexer *p, ghostcell *pgc, double *W, double *WH)
{   
    starttime=pgc->timer();
    
    count=0;
    LOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            W[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*wval[count] + relax4_wg(i,j)*W[IJK];
            WH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*WHval[count] + relax4_wg(i,j)*WH[IJK];
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            W[IJK] = relax4_nb(i,j)*W[IJK];
            WH[IJK] = relax4_nb(i,j)*WH[IJK];
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;		
}

void iowave::P_relax(lexer *p, ghostcell *pgc, double *P)
{
    starttime=pgc->timer();
    FLOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
        {            
            // Zone 2
            if(db<1.0e20)

            P[FIJK] = relax4_nb(i,j)*P[FIJK];
        }
    }	
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::turb_relax_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F)
{
    starttime=pgc->timer();
    
    LOOP
    {
        dg = distgen(p);    
        db = distbeach(p);

		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            F[IJK] = relax4_wg(i,j)*F[IJK];

		}
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            F[IJK] = relax4_nb(i,j)*F[IJK];
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}
