/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_cfd::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    for (int nb = 0; nb < number6DOF; nb++)
        fb_obj[nb]->initialize_cfd(p, a, pgc, pnet);
    
    if(p->Y5>0)
    {
        setup_chrono(p,a,pgc,pnet);
    }
}

void sixdof_cfd::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{}

void sixdof_cfd::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_cfd::setup(lexer *p, fdm *a, ghostcell *pgc)
{
    // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}

void sixdof_cfd::initialize_chrono(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet, int nb)
{
    fb_obj[nb]->ini_fbvel(p,pgc);

    fb_obj[nb]->ray_cast(p,a,pgc);
    // fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Calculate geometrical properties
    fb_obj[nb]->geometry_parameters(p,a,pgc);
    
    // Initialise position of bodies
    fb_obj[nb]->iniPosition_RBM(p,pgc);
    
    // Recalculate distances
    fb_obj[nb]->ray_cast(p,a,pgc);
    // fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Initialise global variables
    fb_obj[nb]->update_fbvel(p,pgc);

    // Initialise floating fields
    ULOOP
    a->fbh1(i,j,k) = fb_obj[nb]->Hsolidface(p,a,1,0,0);

    VLOOP
    a->fbh2(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,1,0);

    WLOOP
    a->fbh3(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,0,1);

    LOOP
    a->fbh4(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,0,0);

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);

    // Print initial body 
    if(p->X50==1)
    fb_obj[nb]->print_vtp(p,pgc);
    
    if(p->X50==2)
    fb_obj[nb]->print_stl(p,pgc);

    // ghostcell update
    pgc->gcdf_update(p,a);
}